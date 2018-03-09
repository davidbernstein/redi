/*
 * Copyright (C) 2004-2018 David Bernstein <david.h.bernstein@gmail.com>
 *
 * This file is part of ReDi.
 *
 * ReDi is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * ReDi is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with ReDi.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "simulation.h"
#include "utility.h"
#include "concentration.h"
#include "parser.h"
#include "constants.h"
#include "realmatrix.h"

#include <fstream>
#include <iostream>
#include <iomanip>

using namespace NAMESPACE;
using namespace std;


void Simulation::Run(const std::string &fileName)
{
    // read input file
    ReadInputFile(fileName);
        
    Initialize();
    OutputEnsembleData();
                
    StopFlag flag;
    
    cout << "number of elements " << mMesh.NumElements() << endl;
    
    PrintSizeOfReactionList();
    //mMesh.PrintVolumeByMaterialNumber();
    //mMesh.PrintSurfaceArea();
    //mMesh.PrintLongestEdgeLength();
    
    //PrintZeroPopulationStop();
    
    mChemicalSystem.PrintReactions();
    
    for (long i = 1; i <= mSimControl.NumberOfRuns(); ++i) {
        flag = Run();
    	
        //PrintDiffusionCount();
        //PrintReactionCount();
        
        OutputEnsembleData();
        
        if (flag == CONTACT_WITH_BOUNDARY) {
        	cout << i << "boundary contact" << endl;
        	return;
        }
    }
    
    
    return;
}



StopFlag Simulation::Run()
{	   	
	++mEnsembleCount;
	Reset();

    // reset random number generator
	mRandomNumberGenerator.Reset();
    
    // initialize output 
    OutputData();
     
    if (mSimControl.StartTime() >= mSimControl.EndTime()) 
    	return END_TIME_REACHED;
    
    // start clock
    mRunStatistics.StartClock();
   
    StopFlag flag = NO_STOP;
	while(flag == NO_STOP) {
		flag = Update();
	}
	
	UpdateRunStatistics();
	
	return flag;
}



StopFlag Simulation::Update()
{	    
	// get reaction
	Reaction *pReaction = mReactionRateTree.GetNextReaction(mRandomNumberGenerator.Random01());
	
	// update current time and test for output and stop
	double deltaT = mReactionRateTree.GetNextReactionTime(mRandomNumberGenerator.Random01());
	if (deltaT > 0.0) {
		OutputData(mCurrentTime + deltaT);
    	if (mCurrentTime + deltaT > mSimControl.EndTime()) {
    		mCurrentTime = mSimControl.EndTime();
    		return END_TIME_REACHED;
    	}
	}
	else {
		return NO_REACTION_POSSIBLE;
	}
	
	// update current time
	mCurrentTime += deltaT;

    // check for stop on contact with boundary
    if (mSimControl.StopOnContactWithBoundary()) {
        if (pReaction->GetElement()->OnBoundary()) 
            return CONTACT_WITH_BOUNDARY;
    }
     
    // perform reaction
    PerformReaction(pReaction);

	// update all affected rates
	UpdateDependentReactionRates(pReaction);
	
	// update tree
	mReactionRateTree.UpdateDependentRates(pReaction);
	
	++mStepCount;
	       
	// check for zero population stop
	if (mSimControl.ZeroPopulationStop()) {
		for (short i = 0; i < mSpeciesName.Size(); ++i) {
			if (mSimControl.ZeroPopulationStop(i) && (mPopulation[i] == 0))
				return ZERO_POPULATION;
		}
	}
	
	
    return NO_STOP;
}



void Simulation::PerformReaction(const Reaction *pReaction)
{	
	if (pReaction->NeighborIndex() < 0)
		PerformChemicalReaction(pReaction);
	else
		PerformDiffusionReaction(pReaction);
	
	return;
}



void Simulation::PerformChemicalReaction(const Reaction *pReaction)
{		
	Element *pE = pReaction->GetElement();
	Occupancy *pOcc = pE->GetOccupancyData();
		
	short index = pReaction->ReactionIndex();
	
	short materialNumber = pE->MaterialNumber();
	
	for (short i = 0; i < mChemicalSystem.NumSpecies(); ++i) {	
		short change = mChemicalSystem.ProductCoefficient(index, i) 
                     - mChemicalSystem.ReactantCoefficient(index, i);
        
		pOcc->IncrementNumOccupants(i, change);
		
		mPopulation[i] += change;
		mMoleculeCount(materialNumber, i) += change;
	}
	
	++mReactionCount[index];
	
	return;
}



void Simulation::PerformDiffusionReaction(const Reaction *pReaction)
{	
	short speciesNumber = pReaction->ReactionIndex();
	short neighborNumber = pReaction->NeighborIndex();
	
	Element *pE = pReaction->GetElement();
	Element *pNeighbor = pE->PN(neighborNumber);
	
	pE->GetOccupancyData()->IncrementNumOccupants(speciesNumber, -1);
	pNeighbor->GetOccupancyData()->IncrementNumOccupants(speciesNumber, 1);
	
	mMoleculeCount(pE->MaterialNumber(), speciesNumber) -= 1;
    mMoleculeCount(pNeighbor->MaterialNumber(), speciesNumber) += 1;
    
    ++mDiffusionCount[speciesNumber];
    
	return;
}



void Simulation::UpdateDependentReactionRates(Reaction *pReaction)
{
	// update dependent reactions
	Array<Reaction*>& graph = pReaction->GetDependencyGraph();
	
	for (short i = 0; i < graph.Size(); ++i) 
		SetReactionRate(graph[i]);
	
	return;
}



void Simulation::UpdateRunStatistics()
{
	mRunStatistics.StopClockAndPrintTime();
	mRunStatistics.AddToTotalNumberOfSteps(mStepCount);
	
	return;
}



void Simulation::Initialize()
{
	// set element data
	if (mSimControl.NeedRefinementData())
		mMesh.SetElementDataSize(4);
	else
		mMesh.SetElementDataSize(3);
	
	if (mSimControl.NeedRefinementData()) {
		mMesh.AddElementData(BISECTION_DATA);
		mMesh.SetBisectionData();
	}
	
	// perform initial mesh refinement, if any
	mMesh.RefineAll(mSimControl.InitialRefinementLevel());
	if (mSimControl.RefinementTestOn())
		mMesh.TestRefinement();
	
	mMesh.AddElementData(OCCUPANCY, NumSpecies());
	
	mMesh.AddElementData(DIFFUSION_FACTOR);
	
	if (mMesh.DefaultElementType() == POLYGON) {
		mMesh.AddElementData(POINT_SET);
		mMesh.ReadElementPointSetData("/Users/dave/Projects/ReDi/p400.dat");
	}
		
   	// set time
    mCurrentTime = mSimControl.StartTime();
    
    // set time for first output of data
    mSimControl.ResetNextOutputTimes();
    
	// set occupancy data
	SetElementOccupancy();
	
	// set diffusion factors
	SetElementDiffusionFactors();
	
	// set reaction factors
	mChemicalSystem.ComputeRateFactors();
		
	// initialize list, tree, and graph
	InitializeReactionList();
	mReactionRateTree.Make(mReactionList);
	InitializeElementReactionSet();
	InitializeDependencyGraph();
	
	// initialize run statistics
	mRunStatistics.Initialize();
	
	// ensemble count
	mEnsembleCount = 0;
	
	
	return;
}



void Simulation::Reset()
{
	// initialize run statistics
	mRunStatistics.Initialize();
	
	// set time
    mCurrentTime = mSimControl.StartTime();
    
    // step count
    mStepCount = 0;
    
    // reaction count
    for (short i = 0; i < mReactionCount.Size(); ++i)
    	mReactionCount[i] = 0;
    	
    // diffusion count
    for (short i = 0; i < mDiffusionCount.Size(); ++i)
    	mDiffusionCount[i] = 0;
    	
    // set time for first output of data
    mSimControl.ResetNextOutputTimes();
    
	// set occupancy data
	SetElementOccupancy();
	
	SetReactionRates();

	mReactionRateTree.SetValues();
	
	InitializeMoleculeCount();

	
	return;
}



void Simulation::InitializeReactionList()
{
	// erase old list
	mReactionList.erase(mReactionList.begin(), mReactionList.end());
	
	Reaction reaction;

	Array<Element*> pE;
	mMesh.GetElements(pE);
    
    short numReactions = mChemicalSystem.NumReactions();
    short numSpecies = mChemicalSystem.NumSpecies();
    
    double startTime = mSimControl.StartTime();
    
    for (long i = 0; i < pE.Size(); ++i) { 
    	short m = pE[i]->MaterialNumber();
    	
    	// add chemical reactions
    	for (short j = 0; j < numReactions; ++j) {
    		if (mChemicalSystem.ReactionActive(j, m)) {
    			reaction.SetReactionIndex(j);
    			reaction.SetNeighborIndex(-1);
    			reaction.SetElement(pE[i]);
    			AddReaction(reaction);
    		}
    	}

    	// add diffusion reactions
    	for (short j = 0; j < numSpecies; ++j) {  
    		if (mMaterial[m].DiffusionCoefficient(j) > 0.0) {
    			for (short k = 0; k < pE[i]->NumNeighbors(); ++k) {
    				if (pE[i]->PN(k)->Type() != BLOB) {
    		   			reaction.SetReactionIndex(j);
    		   			reaction.SetNeighborIndex(k);
    					reaction.SetElement(pE[i]);
    					AddReaction(reaction);
    				}
    			}
    		}
    	}
    }
    
    // set rates
    SetReactionRates();
    
    
	return;
}



void Simulation::SetReactionRates()
{
	list<Reaction>::iterator iRl;
    for (iRl = mReactionList.begin(); iRl != mReactionList.end(); ++iRl)
    	SetReactionRate(&*iRl);
    	
	return;
}



void Simulation::InitializeElementReactionSet()
{
	mElementReactionSet.clear();
	
	list<Reaction>::iterator iRl;
	for (iRl = mReactionList.begin(); iRl != mReactionList.end(); ++iRl) {
		Element *pE = (*iRl).GetElement();
		Reaction *pR = &*iRl;
		mElementReactionSet.insert(make_pair(pE, pR));
	}
	
	return;
}



void Simulation::InitializeDependencyGraph()
{
	list<Reaction>::iterator iRl;
	for (iRl = mReactionList.begin(); iRl != mReactionList.end(); ++iRl) 
		InitializeDependencyGraph(&*iRl);
	
	//CheckDependencyGraph();
	
	return;
}



void Simulation::InitializeDependencyGraph(Reaction *pReaction) 
{
	set<Reaction*> reactionSet;
	
	Element *pE = pReaction->GetElement();
	GetElementReactions(pE, reactionSet, true);
		
	for (short i = 0; i < pE->NumNeighbors(); ++i) 
		GetElementReactions(pE->PN(i), reactionSet, true);
		
	set<Reaction*>::iterator iRs = reactionSet.begin();
	while (iRs != reactionSet.end()) {
		set<Reaction*>::iterator iRsPrev = iRs;
		++iRs;
		if (ReactionADependsOnReactionB(*iRsPrev, pReaction) == false) {
			reactionSet.erase(iRsPrev);
		}
	}

	pReaction->InitializeDependencyGraph(reactionSet);
		
		
	return;
}



bool Simulation::ReactionADependsOnReactionB(const Reaction *pA, const Reaction *pB) const 
{
	short rA = pA->ReactionIndex();
	short nA = pA->NeighborIndex();
	
	short rB = pB->ReactionIndex();
	short nB = pB->NeighborIndex();
	
	Element *pEA = pA->GetElement();
	Element *pEB = pB->GetElement();
	
	if (nA < 0) {
		if (nB < 0) {
			// A is a chemical reaction and B is a chemical reaction
			if (pEA != pEB)
				return false;
			else
				return mChemicalSystem.ReactionDependentOnReaction(rA, rB);
		}
		else {
			// A is a chemical reaction and B is a diffusion reaction
			if (pEA == pEB) {
				return mChemicalSystem.ReactionDependentOnSpecies(rA, rB);
			}
			else {
				if (pEB->NeighborIndex(pEA) != nB)
					return false;
				else
					return mChemicalSystem.ReactionDependentOnSpecies(rA, rB);
			}
		}
	}
	else {
		if (nB < 0) {
			// A is a diffusion reaction and B is a chemical reaction
			if (pEA == pEB) 
				return (mChemicalSystem.ReactantCoefficient(rB, rA) > 0) || (mChemicalSystem.ProductCoefficient(rB, rA) > 0);
			else 
				return false;
		}
		else {
			// A is a diffusion reaction and B is a diffusion reaction
			if (pEA == pEB) {
				return rA == rB;
			}
			else {
				if (pEB->NeighborIndex(pEA) != nB)
					return false;
				else
					return rA == rB;
			}
		}
	}
}



void Simulation::GetElementReactions(const Element *pE, set<Reaction*> &reactionSet, bool addToSet) const
{
	pair<Element*, Reaction*> tmp((Element*) pE, NULL);
	
	ElementReactionSet::iterator iS = mElementReactionSet.lower_bound(tmp);
	
	Element *pEfirst = (*iS).first;
	
	if (addToSet == false)
		reactionSet.clear();
	
	while ((pEfirst == pE) && (iS != mElementReactionSet.end())) {
		reactionSet.insert((*iS).second);
		pEfirst = (*iS).first;
		++iS;
	}
	

	return;
}



void Simulation::SetElementOccupancy()
{
	Array<Element*> pE;
	mMesh.GetElements(pE);

	// erase current occupancy
	for (long j = 0; j < NumSpecies(); ++j) {
		for (long i = 0; i < pE.Size(); ++i) 
			pE[i]->GetOccupancyData()->SetNumOccupants(j, 0);
	}

	for (short j = 0; j < mInitialConcentration.Size(); ++j) {
		Concentration concentration = mInitialConcentration[j];
		
		if (concentration.Type() == CONVERGENCE_TEST)
			concentration.SetParameter(4, mSimControl.StartTime());
		
		if (concentration.Type() == DELTA_FUNCTION) {
			Element *pElement = (Element*) mInitialConcentration.DeltaFunctionElement(j);
		
			long numMolecules = NearestInteger(concentration.Parameter(3));
        	if (numMolecules > MACHINE_LARGEST_SHORT)
            	ThrowException("QMesh::SetElementOccupancy : number of molecules has overflowed short");
        
			if (pElement == NULL) {
				long i = 0;
				bool found = false;
				while ((!found) && (i < pE.Size())) {
            		Point p(concentration.Parameter(0), concentration.Parameter(1), concentration.Parameter(2));
            		found = pE[i]->ContainsPoint(p);
            		if (found) {
						pElement = pE[i];
						mInitialConcentration.SetDeltaFunctionElement(j, pElement);    		
            		}
            		++i;
    			}
    			if (!found)
    				ThrowException("Simulation::SetElementOccupancy : delta function not in domain");
    		}
    		
    		pElement->GetOccupancyData()->SetNumOccupants(j, numMolecules);
		}
		else {
			vector<double> bin(pE.Size());
	
			Point centroid;
			double sum = 0.0;
			for (long i = 0; i < pE.Size(); ++i) {
				pE[i]->Centroid(centroid);
				bin[i] = ConcentrationIntegral(concentration, pE[i], pE[i]->MaterialNumber());		
        		sum += bin[i];
			}
	
			mRandomNumberGenerator.MakeBins(bin);
			
			long totalNumber = mRandomNumberGenerator.FloorOrCeiling(sum);
			for (long i = 0; i < totalNumber; ++i) {
				long k = mRandomNumberGenerator.RandomBin(bin);
				pE[k]->GetOccupancyData()->IncrementNumOccupants(j);
			}
		}
	}
	
	
	return;
}



double Simulation::ConcentrationIntegral(const Concentration &concentration, 
										 const Element *pE, short materialNumber) const
{
	double number = 0.0;
	
	short numVertices = pE->NumVertices();
	
	if (pE->Type() == POLYGON) {
    	Vertex centroid;
    	//pE->Centroid(centroid);
    	
    	PointSet *pPS = pE->GetPointSet();
    	centroid = *(pPS->GetPoint(0));
    		
    	Triangle triangle;
    	for (short j = 0; j < numVertices; ++j) {
    		short k = (j + 1) % numVertices;
    		double sum = concentration.Value(pE->PV(j), materialNumber);
    		sum += concentration.Value(pE->PV(k), materialNumber);
    		sum += concentration.Value(&centroid, materialNumber);
    			
    		triangle.SetVertices(&centroid, pE->PV(j), pE->PV(k));
    		sum *= triangle.Area() / 3.0;
    			
    		number += sum;
    	}
    }
    else {
    	for (short j = 0; j < numVertices; ++j) 
    		number += concentration.Value(pE->PV(j), materialNumber);
    	
    	number *= pE->Volume() / numVertices;
    }
    
	return number;
}



void Simulation::SetElementDiffusionFactors()
{
	Array<Element*> pE;
	mMesh.GetElements(pE);
	
	DiffusionFactor *pDF;
	Array<double> faceArea;
	
	for (long i = 0; i < pE.Size(); ++i) {
		pDF = pE[i]->GetDiffusionFactor();
		pDF->SetSize(pE[i]->NumNeighbors());
    	
    	//Array<Point> areaNormal;
    	//pE[i]->OutwardAreaNormal(areaNormal);
    	
    	pE[i]->FaceArea(faceArea);
    	
    	double volume = pE[i]->Volume();
    	if (volume <= 0.0)
    		ThrowException("Simulation::SetElementDiffusionFactors : non-positive volume");

    	Point c, cn; 
    	if (pE[i]->Type() != POLYGON) {
       		pE[i]->Centroid(c);
    	}
    	else {
    		PointSet *pPS = pE[i]->GetPointSet();
    		c = *(pPS->GetPoint(0));
    	}
    	
		for (short j = 0; j < pE[i]->NumNeighbors(); ++j) {
			if (pE[i]->PN(j)->Type() != BLOB) {
				
				if (pE[i]->PN(j)->Type() != POLYGON) {
       				pE[i]->PN(j)->Centroid(cn);
    			}
    			else {
    				PointSet *pPS = pE[i]->PN(j)->GetPointSet();
    				cn = *(pPS->GetPoint(0));
    			}
    			
    			/*
    			Point q = c - cn;
    			Vertex *pV[2];
    			pE[i]->GetEdge(j, pV);
    			Point qq = (*pV[0]) - (*pV[1]);
    			double xx = q.Dot(qq) / (q.Norm() * qq.Norm());
    			if (abs(xx) > 0.0001)
    				cout << j << " XX " << q.Dot(qq) / (q.Norm() * qq.Norm()) << endl;
    			*/
    			
				pDF->SetFactor(j, faceArea[j] / (volume * c.DistanceToPoint(cn)));
				
				/*
				double dotProduct = areaNormal[j].Dot(cn - c);
				if (dotProduct <= 0.0) {
					areaNormal[j].Print();
					c.Print();
					cn.Print();
					pE[i]->Print();
					ThrowException("Simulation::SetElementDiffusionFactors : non-positive diffusion factor");
				}
 				pDF->SetFactor(j, dotProduct / (volume * c.SqrDistanceToPoint(&cn)));
 				*/
 			}
 			else {
 				pDF->SetFactor(j, -1.0);
 			}
 		}
	}
	
	
	return;
}



void Simulation::OutputData(double time)
{
	double outputTime;
	
	if (mSimControl.Output(time, outputTime, MESH_DATA_OUTPUT)) 
		OutputMeshData(outputTime);
    
    if (mSimControl.Output(time, outputTime, CONSOLE_DATA_OUTPUT)) 
    	OutputConsoleData(outputTime);
    
    if (mSimControl.Output(time, outputTime, TIME_DATA_OUTPUT))   		
    	OutputTimeData(outputTime);
    
    if (mSimControl.Output(mCurrentTime, outputTime, TIME_AVERAGE))   
    	OutputTimeAverageData(outputTime, false);
    	
    if (mSimControl.Output(time, outputTime, TIME_AVERAGE_OUTPUT)) 
    	OutputTimeAverageData(outputTime, true);
    
    if (mSimControl.Output(time, outputTime, ENSEMBLE_AVERAGE_TIME_OUTPUT)) 
    	CollectEnsembleTimeAverageData(outputTime);
    
	return;
}



void Simulation::OutputTimeData(double outputTime)
{
	static ofstream timeFile;

	if (mStepCount == 0) {
    	timeFile.open("/Users/dave/MatlabFiles/BonnerData/timeData.dat");
    }
    
	timeFile << outputTime << " ";
    for (short m = 0; m < mMaterial.Size(); ++m) {
		for (short j = 0; j < NumSpecies(); ++j) {
        	timeFile << mMoleculeCount(m, j) << " ";
    	}
    }
    timeFile << endl;
    	
    if (outputTime >= mSimControl.EndTime()) {
    	timeFile.close();
    }

	return;
}



void Simulation::OutputMeshData(double outputTime)
{
/*
	static XYData xyData(2);
	static XYDisplay xyDisplay1, xyDisplay2;
	static MeshDisplay meshDisplay(DISPLAY_MATERIAL, 2);

	if (mStepCount == 0) { 	
    	xyDisplay1.SetXTitle("x");
    	xyDisplay2.SetXTitle("x");
    	xyDisplay1.SetYTitle("X");
    	xyDisplay2.SetYTitle("Y");
    	
    	xyDisplay1.SetYAxisRange(0, 20.0);
    	xyDisplay2.SetYAxisRange(0, 20.0);
    	
    	xyDisplay1.SetXAxisRange(0, 1.0);
    	xyDisplay2.SetXAxisRange(0, 1.0);
    }
    
	if (mMesh.Dimension() == 1) {
        xyData.Make(mMesh, 2);
		xyDisplay1.Show(xyData);
			
		xyData.Make(mMesh, 3);
		xyDisplay2.Show(xyData);
	}
	else {
		meshDisplay.ShowMesh(mMesh);
	}
*/

	return;
}



void Simulation::OutputTimeAverageData(double outputTime, bool output)
{
	static Matrix<double> aveConc;
    static ofstream aveConcFile;
    static long aveCount;
            
    if (mStepCount == 0) {
    	if (mSimControl.OutputOn(TIME_AVERAGE_OUTPUT)) {
    		aveConc.SetSize(mChemicalSystem.NumSpecies(), mMesh.NumElements());
    		aveConcFile.open("/Users/dave/Projects/ReDi/aveConc.dat");
    		aveConcFile << setprecision(8);
    		aveCount = 0;
    	}
    }
    
	if (output) {
    	for (long j = 0; j < aveConc.NumColumns(); ++j) {		
    		for (short i = 0; i < aveConc.NumRows(); ++i) {
    			aveConcFile << aveConc(i, j) << " ";
    		}
    		aveConcFile << endl;
    	}
    	aveCount = 0;
	}
	else {
    	short neighborhoodSize = mSimControl.SpatialAverageNeighborhoodSize();
    
    	if (neighborhoodSize > 1)
    		mMesh.SetVertexElementNeighbors();
    
		Array<Element*> pE;
		mMesh.GetElements(pE);
		
		double factor = 1.0 / (aveCount + 1.0);
		
		set<Element*> neighborhood;
		set<Element*>::iterator iNs;
		
    	for (short i = 0; i < mChemicalSystem.NumSpecies(); ++i) {
    		for (long j = 0; j < pE.Size(); ++j) {
    			pE[j]->Neighborhood(neighborhood, neighborhoodSize, mMesh.DefaultElementType());
    			double sum = 0.0;
    			
    			for (iNs = neighborhood.begin(); iNs != neighborhood.end(); ++iNs) 
    				sum += (*iNs)->GetOccupancyData()->NumOccupants(i);
    			
    			sum /= neighborhood.size();
    			
    			aveConc(i, j) = aveConc(i, j) * aveCount * factor + sum * factor;
    		}
    	}
    	
    	++aveCount;
    	
    	if (neighborhoodSize > 1)
    		mMesh.EraseVertexElementNeighbors();
    }
    
    if (outputTime >= mSimControl.EndTime()) {
    	aveConcFile.close();
    }
    
	return;
}



void Simulation::OutputConsoleData(double outputTime)
{
	cout << "Time " << outputTime << endl;
    
    PrintStepCount();
    PrintMoleculeCount();
    PrintReactionCount();
    PrintDiffusionCount();

	return;
}



void Simulation::CollectEnsembleTimeAverageData(double outputTime, bool dumpData)
{
	static Array<Matrix<long> > population;
	static Array<double> time;
    static long count;
    short numSpecies, numMaterials;
	
    if (mEnsembleCount == 0) {
    	numSpecies = mChemicalSystem.NumSpecies();
    	numMaterials = mMaterial.Size();
    	
    	double interval = mSimControl.OutputInterval(ENSEMBLE_AVERAGE_TIME_OUTPUT);
    	long s = NearestInteger((mSimControl.EndTime() - mSimControl.StartTime()) / interval) + 1;
    	population.SetSize(s);
    	time.SetSize(s);
    	for (long i = 0; i < s; ++i) {
    		population[i].SetSize(numMaterials, numSpecies);
    	
    		for (short m = 0; m < numMaterials; ++m) {
    			for (short s = 0; s < numSpecies; ++s) {
    				population[i](m, s) = 0;
    			}
    		}
    	}
    	
    	return;
    }
    
     
	if (dumpData) {
    	ofstream avePopFile("/Users/dave/Projects/ReDi/avePop.dat");
    	avePopFile << setprecision(8);
    		
    	for (long i = 0; i < count; ++i) {
    		avePopFile << time[i] << " " << (double) population[i](1, 1) / mEnsembleCount << endl;
    	}
    		
    	avePopFile.close();
    	
    	return;
	}
	
    if (mStepCount == 0) 
    	count = 0;
	
	population[count] += mMoleculeCount;
	time[count] = outputTime;
	++count;
	
	
	return;
}



void Simulation::OutputEnsembleData()
{
	if (mSimControl.OutputEndTime())
		OutputEndTime();
	
	if (mSimControl.ConvergenceTestErrorOn())
		DiffusionTestError();
	
	if (mSimControl.RandomWalkTestOn())
		RandomWalkError();
    
	if (mSimControl.CubicTestErrorOn())
		CubicWaveError();
	
	CollectEnsembleTimeAverageData(-1.0, true);


	return;
}



void Simulation::ReadInputFile(const string &fileName)
{
	Parser parser(fileName);
	
	// number of species
	long numSpecies;
	if (parser.FindInteger("numberofspecies=", numSpecies) == false)
		ThrowException("Simulation::ReadInputFile : file " + fileName + " does not contain number of species");
	SetNumSpecies(numSpecies);
	mDiffusionCount.SetSize(numSpecies);
	
	// number of reactions
	long numReactions = 0;
	parser.FindInteger("numberofreactions=", numReactions);
	mChemicalSystem.SetNumReactionsAndSpecies(numReactions, numSpecies);
	mReactionCount.SetSize(numReactions);
	
	// number of materials
	long numMaterials;
	if (parser.FindInteger("numberofmaterials=", numMaterials) == false)
		ThrowException("Simulation::ReadInputFile : didn't find number of materials");
	SetNumMaterials(numMaterials);
		
	// species names
	string n;
	string prefix = "nameofspecies";
	for (short i = 1; i <= numSpecies; ++i) {
		n = ConvertIntegerToString(i);
		string name;
		if (parser.FindString(prefix + n + "=", name) == false)
			ThrowException("Simulation::ReadInputFile : didn't find species name " + n);
		else
			SetSpeciesName(i - 1, name);
	}
	mChemicalSystem.SetSpeciesNames(mSpeciesName);
	
	// units
	string units;
	if (parser.FindString("lengthunits=", units))
		mUnits.SetLengthUnits(units);
	else
		mUnits.SetLengthUnits("meters");
	
	if (parser.FindString("timeunits=", units))
		mUnits.SetTimeUnits(units);
	else
		mUnits.SetTimeUnits("seconds");
		
	// start and end times
	double startTime, endTime;
	if (parser.FindFloat("starttime=", startTime) == false)
		ThrowException("Simulation::ReadInputFile : file " + fileName + " does not contain start time");

	if (parser.FindFloat("endtime=", endTime) == false)
		ThrowException("Simulation::ReadInputFile : file " + fileName + " does not contain end time");
	
	mSimControl.SetRunTime(startTime, endTime);
	
	// diffusion coefficients
	prefix = "diffusioncoefficient";
	Array<double> parameter(numSpecies);
	for (short i = 0; i < numMaterials; ++i) {
		n = ConvertIntegerToString(i+1);
		if (parser.FindBracedFloats(prefix + n + "={", parameter) == false) {
			ThrowException("Simulation::ReadInputFile : didn't find diffusion coefficients for material " + n);
		}
		else {
			mMaterial[i].SetDiffusionCoefficients(parameter);
		}
	}
	
	// reaction rates
	prefix = "reactionrate";
	for (short i = 1; i <= numReactions; ++i) {
		double rate;
		n = ConvertIntegerToString(i);
		if (parser.FindFloat(prefix + n + "=", rate) == false)
			ThrowException("Simulation::ReadInputFile : didn't find reaction rate " + n);
		else
			mChemicalSystem.SetRateConstant(i - 1, rate);
	}
	
	// reaction coefficients
    parameter.SetSize(NumSpecies());
    Array<short> coef(NumSpecies());
    for (short i = 0; i < numReactions; ++i) {
        n = ConvertIntegerToString(i + 1);
        if (parser.FindBracedFloats("reactantcoefficients" + n + "={", parameter) == false)
            ThrowException("Simulation::ReadInputFile : didn't find all reactant coefficients" + n);
            
        for (short j = 0; j < parameter.Size(); ++j)
            coef[j] = (short) NearestInteger(parameter[j]);
            
        mChemicalSystem.SetReactantCoefficients(i, coef);
        
        if (parser.FindBracedFloats("productcoefficients" + n + "={", parameter) == false)
            ThrowException("Simulation::ReadInputFile : didn't find all product coefficients" + n);
        
        for (short j = 0; j < parameter.Size(); ++j)
            coef[j] = (short) NearestInteger(parameter[j]);
            
        mChemicalSystem.SetProductCoefficients(i, coef);
    }
    
    // reaction material activity
    mChemicalSystem.InitializeReactionMaterialActivity(numMaterials);
    parameter.SetSize(numMaterials);
    Array<bool> activity(numMaterials);
    for (short i = 0; i < numReactions; ++i) {
        n = ConvertIntegerToString(i + 1);
        parser.FindBracedFloats("reactionmaterialactivity" + n + "={", parameter);
            
        for (short j = 0; j < parameter.Size(); ++j)
            activity[j] = NearestInteger(parameter[j]) == 1 ? true : false;
            
        mChemicalSystem.SetReactionMaterialActivity(i, activity);
    }
    
	// initial concentration functions
	mInitialConcentration.SetSize(numSpecies);
	prefix = "concentrationtype";
	for (short i = 0; i < numSpecies; ++i) {
		n = ConvertIntegerToString(i + 1);
		string name;
		if (parser.FindString(prefix + n + "=", name) == false)
			ThrowException("Simulation::ReadInputFile : didn't find concentration " + n);
		else 
			mInitialConcentration[i].SetType(name);
		
		parameter.SetSize(mInitialConcentration[i].NumParameters());
		
		if (parser.FindBracedFloats("concentration" + n + "={", parameter) == false)
			ThrowException("Simulation::ReadInputFile : didn't find all parameters for concentration " + n);
			
		for (short j = 0; j < parameter.Size(); ++j) 
			mInitialConcentration[i].SetParameter(j, parameter[j]);
	}    

	// mesh refinement
	string dum;
	if (parser.FindString("meshrefinement=on", dum))
		mSimControl.TurnOnMeshRefinement();
		
	// mesh file name
	string meshFileName;
	if (parser.FindFileName("mesh=", meshFileName) == false)
		ThrowException("Simulation::ReadInputFile : didn't find mesh file");
	
	mMesh.ReadFromFile(meshFileName);
	
	if (mMesh.CheckMaterialNumbers() > numMaterials)
		ThrowException("Simulation::ReadInputFile : mesh has more materials than simulation file");
	
    // stop on contact with boundary
    if (parser.FindString("stoponcontactwithboundary=on", dum))
		mSimControl.SetToStopOnContactWithBoundary();
    
    // stop on population count
    mSimControl.SetStopOnZeroPopulationSize(numSpecies);
    for (short j = 0; j < numSpecies; ++j) {
    	string name = "stopwhenpopulationzero=" + mSpeciesName[j];
    	LowerCase(name);
    	if (parser.FindString(name, dum)) {
			mSimControl.SetToStopOnZeroPopulation(j);
		}
	}
	
    // number of runs (default is 1)
    long numRuns;
	if (parser.FindInteger("numberofruns=", numRuns))
		mSimControl.SetNumberOfRuns(numRuns);
        
    // initial random seed (default is 1)
    long initialSeed = -1;
	parser.FindInteger("randomnumberseed=", initialSeed);
	mRandomNumberGenerator.Reset(-abs(initialSeed));
     
    string rngType;
    if (parser.FindString("randomnumbergenerator=", rngType))
        ThrowException("code not implemented");
        
    // data write interval
    double interval;
    if (parser.FindFloat("meshoutputinterval=", interval))
		mSimControl.SetOutputInterval(interval, MESH_DATA_OUTPUT);
    
    if (parser.FindFloat("timeoutputinterval=", interval))
		mSimControl.SetOutputInterval(interval, TIME_DATA_OUTPUT);
    
    if (parser.FindFloat("consoleoutputinterval=", interval))
		mSimControl.SetOutputInterval(interval, CONSOLE_DATA_OUTPUT);
    
    if (parser.FindFloat("ensembleoutputinterval=", interval))
		mSimControl.SetOutputInterval(interval, ENSEMBLE_DATA_OUTPUT);
    
    if (parser.FindFloat("ensembleaveragetimeinterval=", interval))
		mSimControl.SetOutputInterval(interval, ENSEMBLE_AVERAGE_TIME_OUTPUT);
    
    // initial mesh refinement
    long initalLevel;
    if (parser.FindInteger("initialrefinementlevel=", initalLevel))
		mSimControl.SetInitialRefinementLevel((short) initalLevel);
    
    if (parser.FindString("refinementtest=on", dum))
		mSimControl.TurnOnRefinementTest();
    
    if (parser.FindString("refinementtest=on", dum))
		mSimControl.TurnOnUnRefinementTest();
  
  	// error output
  	if (parser.FindString("convergencetesterror=on", dum))
		mSimControl.TurnOnConvergenceTestError();
  	
  	if (parser.FindString("cubictesterror=on", dum))
		mSimControl.TurnOnCubicTestError();
  	
  	if (parser.FindString("randomwalktest=on", dum))
		mSimControl.TurnOnRandomWalkTest();
	
  	if (parser.FindString("outputendtime=on", dum))
		mSimControl.TurnOnEndTimeOutput();
  	
  	// time averaging
    if (parser.FindFloat("timeaverageoutputinterval=", interval)) {
		mSimControl.SetOutputInterval(interval, TIME_AVERAGE_OUTPUT);
    
  		long idum;
  		if (parser.FindInteger("numberoftimeaveragesperoutput=", idum) == false) 
			idum = mMesh.NumElements();
			
		mSimControl.SetOutputInterval(interval / idum, TIME_AVERAGE);
  	}
  	
  	// space averaging
  	long idum;
  	if (parser.FindInteger("spatialaverageneighborhoodsize=", idum))
		mSimControl.SetSpatialAverageNeighborhoodSize(idum);
  	
	return;
}



void Simulation::InitializeMoleculeCount() 
{
	mMoleculeCount.SetSize(mMaterial.Size(), NumSpecies());
	
	Array<Element*> pE;
	mMesh.GetElements(pE);
   
   	for (short j = 0; j < NumSpecies(); ++j) {
		for (short m = 0; m < mMaterial.Size(); ++m) {
        	mMoleculeCount(m, j) = 0;
        	for (long i = 0; i < pE.Size(); ++i) {
        		if (pE[i]->MaterialNumber() == m) {
            		mMoleculeCount(m, j) += pE[i]->GetOccupancyData()->NumOccupants(j);
        		}
            }
    	}
    }
    
   	for (short j = 0; j < NumSpecies(); ++j) {
   		mPopulation[j] = 0;
		for (short m = 0; m < mMaterial.Size(); ++m) 
        	mPopulation[j] += mMoleculeCount(m, j);
    }
    
    
	return;
}



void Simulation::CenterOfMass(Array<Point> &centerOfMass)
{
	centerOfMass.SetSize(NumSpecies());
	
	Array<Element*> pE;
    mMesh.GetElements(pE);
   
   	Point c;
   	
	for (short j = 0; j < NumSpecies(); ++j) {
        centerOfMass[j].SetPosition(0.0, 0.0, 0.0);
        long count = 0;
        
        for (long i = 0; i < pE.Size(); ++i) {
        	pE[i]->Centroid(c);
            short n = pE[i]->GetOccupancyData()->NumOccupants(j);
            if (n > 0) {
            	centerOfMass[j] += c * n;
            	count += n;
            }
        }
        
        centerOfMass[j] /= (double) count;
    }
    
	return;
}



void Simulation::CheckDependencyGraph() const
{
	// warning: this is an R^2 operation (R is the size of mReactionList)
	// use with caution
	
	list<Reaction>::const_iterator iRl1, iRl2;
	for (iRl1 = mReactionList.begin(); iRl1 != mReactionList.end(); ++iRl1) {
		for (iRl2 = mReactionList.begin(); iRl2 != mReactionList.end(); ++iRl2) {
			if (ReactionADependsOnReactionB(&*iRl1, &*iRl2)) {
				if ((*iRl2).IsADependentReaction(&*iRl1) == false) {
					(*iRl1).Print();
					(*iRl2).Print();
					ThrowException("Simulation::CheckDependencyGraph : bad graph");
				}
			}
			else {
				if ((*iRl2).IsADependentReaction(&*iRl1) == true) {
					(*iRl1).Print();
					(*iRl2).Print();
					ThrowException("Simulation::CheckDependencyGraph : bad graph");
				}
			}
		}
	}
	
	return;
}



void Simulation::CheckReactionList() const
{
	list<Reaction>::const_iterator iRl;
	for (iRl = mReactionList.begin(); iRl != mReactionList.end(); ++iRl) {
		Element *pE = (*iRl).GetElement();
		
		Occupancy *pOcc = pE->GetOccupancyData();
		
		BisectionData *pBD = pE->GetBisectionData();
	}
	
	return;
}



void Simulation::PrintCurrentTime() const
{
    cout << "Time " << mCurrentTime << endl;
    return;
}



void Simulation::PrintStepCount() const
{
	string s = ConvertIntegerToString(mStepCount, true);
	cout << "Step count " << s << endl;
	
	return;
}



void Simulation::PrintMoleculeCount() const
{
	for (short j = 0; j < NumSpecies(); ++j) {
    	cout << "Total number of " + mSpeciesName[j] + " = "  << mPopulation[j] << endl;
   	}
   	
   	cout << endl;
   	
   	for (short m = 0; m < mMaterial.Size(); ++m) {
   	   string matString = " in material " + ConvertIntegerToString(m);
       for (short j = 0; j < NumSpecies(); ++j) {
       	   cout << "Number of " + mSpeciesName[j] + matString + " = " << mMoleculeCount(m, j) << endl;
   	   }
    }
    
    return;
}



void Simulation::PrintReactionCount() const
{
	for (short i = 0; i < mReactionCount.Size(); ++i) {
		cout << "Reaction count " << i << " = ";
		cout << ConvertIntegerToString(mReactionCount[i], true);
		cout << endl;
	}
	
	cout << endl;

	return;
}



void Simulation::PrintDiffusionCount() const
{
	for (short i = 0; i < mDiffusionCount.Size(); ++i) {
		cout << "Diffusion count " << i << " = ";
		cout << ConvertIntegerToString(mDiffusionCount[i], true);
		cout << endl;
	}
	
	cout << endl;

	return;
}



void Simulation::PrintSizeOfReactionList() const
{
	cout << "Size of reaction list " << mReactionList.size() << endl;

	return;
}



void Simulation::CubicWaveError() const
{
	static Array<double> averageError, averageSolution;
	static double outputTime, firstOutputTime = 1.0;
	static ofstream file, r2file;
	
	Array<Element*> pE;
    mMesh.GetElements(pE);
    long numElements = pE.Size();
    
    if (mEnsembleCount == 0) {
    	averageError.SetSize(numElements);
    	averageSolution.SetSize(numElements);
    	
    	outputTime = firstOutputTime;
    	
    	file.open("/Users/dave/Projects/ReDi/cubicerror.dat");
    	
    	for (long i = 0; i < numElements; ++i) 
    		averageError[i] = 0.0; 
    		
    	return;
    }
    
    double factor = ((double) mEnsembleCount - 1.0) / (mEnsembleCount);
       
   	double maxError = -1.0, aveError = 0.0;
   	   	   	
   	Concentration concentration = mInitialConcentration[1];
   	concentration.SetParameter(3, mCurrentTime);
   	
   	long iMax = 0;
   	
   	double gridFraction = 0.617;
   	long numE = NearestInteger(gridFraction * pE.Size());
   	
    for (long i = 0; i < numE; ++i) {
    	double exact = 0.0;
    	for (short j = 0; j < pE[i]->NumVertices(); ++j) 
    		exact += concentration.Value(pE[i]->PV(j));
    	
    	exact *= pE[i]->Volume() / pE[i]->NumVertices();
    			
		double sol = pE[i]->GetOccupancyData()->NumOccupants(1);
				
		double error = 1.0 - sol / exact;
		
    	averageError[i] = averageError[i] * factor + error / mEnsembleCount;
    	
		aveError += fabs(averageError[i]);
		
		if (fabs(averageError[i]) > maxError) {
			maxError = fabs(averageError[i]);
			iMax = i;
		}
    }
    aveError /= numE;
    
    for (long i = 0; i < numElements; ++i) {
    	double sol = pE[i]->GetOccupancyData()->NumOccupants(1) / pE[i]->Volume();
		averageSolution[i] = averageSolution[i] * factor + sol / mEnsembleCount; 
    }
    
    cout << mEnsembleCount << " " << aveError << " " << maxError << " " << iMax << endl;
    
    double logCount = log10((double) mEnsembleCount);
    if (logCount >= firstOutputTime) {
    	if ((logCount > outputTime) || (mEnsembleCount == mSimControl.NumberOfRuns())) {
    		file << logCount << " " << log10(fabs(aveError)) <<  " ";
    		file << log10(maxError) << endl;
    	
    		outputTime += mSimControl.OutputInterval(ENSEMBLE_DATA_OUTPUT);
    		
    		cout << mEnsembleCount << " " << aveError << " " << maxError << endl;
    	}
    }
    
    if (mEnsembleCount == mSimControl.NumberOfRuns()) {
    	Point centroid;
    	ofstream aveFile("/Users/dave/Projects/ReDi/avecubicsol.dat");
    	double a0 = mInitialConcentration[1].Parameter(0);
    	for (long i = 0; i < numElements; ++i) {
    		pE[i]->Centroid(centroid);
    		aveFile << centroid.X() << " " << averageSolution[i] / a0 << endl;
    	}
    	
    	ofstream typFile("/Users/dave/Projects/ReDi/typicalsol.dat");
    	for (long i = 0; i < numElements; ++i) {
    		pE[i]->Centroid(centroid);
    		double sol = pE[i]->GetOccupancyData()->NumOccupants(1) / pE[i]->Volume();
    		typFile << centroid.X() << " " << sol / a0 << endl;
    	}
    }
    
	return;
}



void Simulation::DiffusionTestError() const
{
	static double maxError, outputTime, outputInterval = 0.1;
	static Array<double> averageError, averageSolution;
    static ofstream file;
	
	Array<Element*> pE;
    mMesh.GetElements(pE);
    
	if (mEnsembleCount == 0) {
    	averageError = 0.0; 
    	outputTime = 1.0;
    	
    	file.open("/Users/dave/Projects/ReDi/differror.dat");
    	
    	averageError.SetSize(pE.Size());
    	averageSolution.SetSize(pE.Size());
    	for (long i = 0; i < pE.Size(); ++i) {
    		averageError[i] = 0.0; 
    		averageSolution[i] = 0.0;
    	}
    		
    	return;
    }
    
    double factor = ((double) mEnsembleCount - 1.0) / mEnsembleCount;
   	 	
   	double aveErr = 0.0;
   	double maxErr = -1.0;
   	long iMax = 0;
   	
   	Concentration concentration = mInitialConcentration[0];
   	
   	switch (concentration.Type()) {
   	case CONVERGENCE_TEST:
   		concentration.SetParameter(4, mCurrentTime);
   		break;
   	
   	case DISCONTINUITY_TEST:
   		concentration.SetParameter(1, mCurrentTime);
   		break;
   		
   	default:
   		ThrowException("Simulation::DiffusionTestError : bad concentration type");
   		break;
   	}
   	
    for (long i = 0; i < pE.Size(); ++i) {
    	double exact = ConcentrationIntegral(concentration, pE[i]);	
		double sol = pE[i]->GetOccupancyData()->NumOccupants(0);
		double error = 1.0 - sol / exact;
		
    	averageError[i] = averageError[i] * factor + error / mEnsembleCount;
    	
		aveErr += fabs(averageError[i]);
		
		if (fabs(averageError[i]) > maxErr) 
			iMax = i;
		
		maxErr = max(fabs(averageError[i]), maxErr);
		
		averageSolution[i] = averageSolution[i] * factor + (sol / pE[i]->Volume()) / mEnsembleCount; 
    }
        
    aveErr /= pE.Size();
    
    cout << mEnsembleCount << " " << log10(fabs(aveErr)) << " " << log10(fabs(maxErr)) << endl;

    double logCount = log10((double) mEnsembleCount);
    if ((logCount >= outputTime) || (mEnsembleCount == mSimControl.NumberOfRuns())) {
    	file << logCount << " " << log10(fabs(aveErr)) <<  " ";
    	file << log10(fabs(maxErr)) << endl;
    	    	
    	outputTime += outputInterval;
    }
    
    if (mEnsembleCount == mSimControl.NumberOfRuns()) {
    	Point centroid;
    	ofstream aveFile("/Users/dave/Projects/ReDi/avediffsol.dat");
    	//double a0 = mInitialConcentration[1].Parameter(0);
    	for (long i = 0; i < pE.Size(); ++i) {
    		pE[i]->Centroid(centroid);
    		aveFile << centroid.X() << " " << averageSolution[i] << endl;
    	}
    	/*
    	ofstream typFile("/Users/dave/Projects/ReDi/typicalsol.dat");
    	for (long i = 0; i < numElements; ++i) {
    		pE[i]->Centroid(centroid);
    		double sol = pE[i]->GetOccupancyData()->NumOccupants(1) / pE[i]->Volume();
    		typFile << centroid.X() << " " << sol / a0 << endl;
    	}
    	*/
    }
    
    
	return;
}



void Simulation::RandomWalkError() const
{
	static double averageDisp, outputTime, outputInterval, aveDimension;
    static ofstream file, r2file;
    if (mEnsembleCount == 0) {
    	averageDisp = 0.0; 
    	aveDimension = 0.0;
    	outputTime = 0.0;
    	outputInterval = 0.2;
    	
    	file.open("/Users/dave/Projects/ReDi/rwmean.dat");
    	r2file.open("/Users/dave/Projects/ReDi/r2.dat");
    	
    	return;
    }
    
	Array<Element*> pE;
	mMesh.GetElements(pE);
	
	long i = 0;
	Element *pE0 = NULL;
	while((pE0 == NULL) && (i < pE.Size())) {
		if (pE[i]->GetOccupancyData()->NumOccupants(0) != 0)
			pE0 = pE[i];
		
		++i;
	}
	
	if (pE0 == NULL)
		ThrowException("Simulation::RandomWalkError : didn't find element");
	
	// make sure element is occupied
    if (pE0->GetOccupancyData()->NumOccupants(0) != 1)
		ThrowException("Simulation::RandomWalkError : bad element occupancy");
			
   	const Element *pEStart = mInitialConcentration.DeltaFunctionElement(0);
   	Point c0, c1;
   	pEStart->Centroid(c0);
   	pE0->Centroid(c1);
   	double r2 = c0.SqrDistanceToPoint(&c1);

   	double bigD = mMaterial[0].DiffusionCoefficient(0);

    double disp = r2 / (2 * bigD * pE0->Dimension() * mCurrentTime);
    double d = r2 / (2 * bigD * mCurrentTime);
    
    averageDisp = averageDisp * (mEnsembleCount - 1.0) / mEnsembleCount + disp / mEnsembleCount; 
    //aveDimension =  aveDimension * n / (n + 1.0) + d / (n + 1.0); 
    
    r2file << r2 << endl;
    
    double logCount = log10((double) mEnsembleCount);
    if (logCount >= 1.0) {
    	if ((logCount > outputTime) || (mEnsembleCount == mSimControl.NumberOfRuns())) {
    		file << logCount << " " << averageDisp <<  " ";
    		file << aveDimension << endl;
    	
    		outputTime += outputInterval;
    		
    		cout << mEnsembleCount << " " << aveDimension << endl;
    	}
    }
        
	//cout << count << " " << error << " " << averageError << endl;
	//c1.Print();
	
	return;
}



void Simulation::OutputEndTime()
{
	static ofstream endTimeFile;
  
    if (mEnsembleCount == 0) {
    	endTimeFile.open("/Users/dave/Projects/ReDi/endTime.dat");
    	endTimeFile << setprecision(12);
    	
    	return;
    }
    
	endTimeFile << (mCurrentTime - mSimControl.EndTime()) / mSimControl.EndTime() << endl;
	    
	return;
}



void Simulation::PrintZeroPopulationStop() const
{
	cout << endl;
	cout << "Zero population stop on = " << mSimControl.ZeroPopulationStop() << endl;
	for (short i = 0; i < mSpeciesName.Size(); ++i) {
		 cout << "species index = " << i << ": " << mSimControl.ZeroPopulationStop(i) << endl;
	}
	
	cout << endl;
	
	return;
}



void Simulation::Test(const string &fileName)
{
	// read input file
    ReadInputFile(fileName);
        
    Initialize();
	
	long numElements = mMesh.NumElements();
	long numReactions = mReactionList.size();
	short treeDepth = mReactionRateTree.Depth();
	cout << "num elements = " << numElements << endl;
	cout << "num reactions = " << numReactions << endl;
	cout << "tree depth = " << treeDepth << endl;
	PrintMoleculeCount();
	
	long n = 10000000;
	long m = 5;
	
	Clock clock;
	double urn1 = 0.0, urn2 = 0.0, urn3 = 0.0, x = 0.0, value = 0.0;
	double y = 1;
	///*
	for (short depth = 0; depth <= treeDepth; ++depth) {
	double sum = 0.0;
	for (short j = 0; j < m; ++j) {
		clock.Start();
		
		long leftCount = 0, rightCount = 0;
		
		for (long i = 0; i < n; ++i) {
			urn1 = 1.0 * mRandomNumberGenerator.Random01();
			//urn2 = 1.0 - ((double) i) / n;
			
			//cout << urn1 << " " << urn2 << endl;
			//Reaction *pReaction = mReactionRateTree.GetNextReaction(urn1);
			//y += urn1;
			///*
			/*
			short d = 0;
			ReactionRateNode* pNode = mReactionRateTree.Root();
			x = urn1 * pNode->Value();
			while ((pNode->IsALeaf() == false) && (d < depth)) {
				value = pNode->LeftChild()->Value();
				
				if (x < value) {
					pNode = pNode->LeftChild();
					//++leftCount;
				}
				else {
					x -= value;
					pNode = pNode->RightChild();
					//++rightCount;
				}
				
				++d;
			}
			*/
			//*/
		}
		
		//cout << "left = " << leftCount << ", right = " << rightCount << " " << (double) leftCount / rightCount << endl;
		sum += clock.StopAndPrintTime("time ");
	}
	//*/
	/*
	list<Reaction>::iterator iRl;
	for (short j = 0; j < m; ++j) {
		StartClock();
		short x = 0;
		for (iRl = mReactionList.begin(); iRl != mReactionList.end(); ++iRl) {
			//SetReactionRate(&*iRl);
			//UpdateDependentReactionRates(&*iRl);
			//x = max(x, (*iRl).NumDependentReactions());
			mReactionRateTree.UpdateDependentRates(&*iRl);
		}
		
		sum += StopClockAndPrintTime("time ");
		cout << x<< endl;
		//sum += StopClock();
	}
	*/
	
	cout << "depth = " << depth << endl;
	double averageTime = sum / (m);
	cout << "average time = " << averageTime  << endl;
	//cout << averageTime / log10((double) numElements) << endl;
	cout << (averageTime - 0.74) / depth << endl;
	}
	
	return;
}
