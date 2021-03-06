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

#ifndef _simulation_h_
#define _simulation_h_

#include "material.h"
#include "reditypes.h"
#include "simcontrol.h"
#include "units.h"
#include "concentrationset.h"
#include "randomnumbergenerator.h"
#include "reaction.h"
#include "multilevelmesh.h"
#include "reactionratetree.h"
#include "runstatistics.h"

namespace NAMESPACE {
	class Simulation {
    public:
        Simulation(void) { };
		~Simulation(void) { };
        
        // run
        void Run(const std::string &fileName);
        
        void Test(const std::string &fileName);
        
	private:
		// run
        StopFlag Run(void);
        StopFlag Update(void);
		void Reset(void);   
		
		// reaction
		void PerformChemicalReaction(const Reaction *pReaction);
		void PerformReaction(const Reaction *pReaction);
		void SetReactionRate(Reaction *pReaction) const;
		
		// diffusion
		void PerformDiffusionReaction(const Reaction *pReaction);
		double DiffusionRate(const Element* pE, short speciesNumber, short neighborNumber) const;
		void SetElementDiffusionFactors(void);
		
		// Reaction list
		Reaction* AddReaction(const Reaction &reaction);
		void InitializeReactionList(void);
		
		// Dependency graph
		void UpdateDependentReactionRates(Reaction *pReaction);
		void InitializeDependencyGraph(void);
		void InitializeDependencyGraph(Reaction *pReaction);
		bool ReactionADependsOnReactionB(const Reaction *pA, const Reaction *pB) const;

		// element-reaction set
		void InitializeElementReactionSet(void);
		void GetElementReactions(const Element *pE, std::set<Reaction*> &reactionSet, bool addToSet = false) const;

		// initialization
		void ReadInputFile(const std::string &fileName);
		void Initialize(void);
        void SetElementOccupancy(void);
        void SetReactionRates(void);
        
        // run statistics
        void UpdateRunStatistics(void);
        
        // species names
        void SetSpeciesName(short index, const std::string &name);
     	std::string SpeciesName(short index) const;
        void SetNumSpecies(short num);
        short NumSpecies(void) const;
        
        // materials
        void SetNumMaterials(short number);
       
        // checking
        void CheckDependencyGraph(void) const;
        void CheckReactionList(void) const;
        
    private:
    	// statistics
		void InitializeMoleculeCount(void);  
		void CenterOfMass(Array<Point> &centerOfMass);
		void DiffusionTestError(void) const;
		void RandomWalkError(void) const;
		void CubicWaveError(void) const;
	 
        // IO
        void OutputData(double time = 0.0);
        void OutputConsoleData(double outputTime);
        void OutputMeshData(double outputTime);
        void OutputTimeData(double outputTime);
        void OutputTimeAverageData(double outputTime, bool output);
        void CollectEnsembleTimeAverageData(double outputTime, bool dumpData = false);
        void OutputEnsembleData(void);
        void OutputEndTime(void);
        void PrintMoleculeCount(void) const;
		void PrintStepCount(void) const;
		void PrintCurrentTime(void) const;
		void PrintReactionCount(void) const;
		void PrintDiffusionCount(void) const;
		void PrintZeroPopulationStop(void) const;
		void PrintSizeOfReactionList(void) const;
		
		// concentration functions
		double ConcentrationIntegral(const Concentration &concentration,
									 const Element *pE, short materialNumber = 0) const;
		
	private:
		// material properties
		Array<Material> mMaterial;
		
		// chemical reactions
		ChemicalSystem mChemicalSystem;
		
        // run controller
        SimControl mSimControl;
        
        // random number generator
        RandomNumberGenerator mRandomNumberGenerator;

		// initial concentration functions
		ConcentrationSet mInitialConcentration;
		
		// list of reactions
		std::list<Reaction> mReactionList;
		
		// reaction rate tree
		ReactionRateTree mReactionRateTree;
		
		// mesh
		MultiLevelMesh mMesh;
		
		// element-reaction set
		ElementReactionSet mElementReactionSet;
		
        // units
        Units mUnits;
        
        // time of run
        double mCurrentTime;
        
        // step count
        long long mStepCount;
        
        // ensemble count
        long mEnsembleCount;
        
        // species names
        Array<std::string> mSpeciesName;
        
        //  reaction count
        Array<long> mReactionCount;
        
        // diffusion count
        Array<long> mDiffusionCount;
        
        // molecule count
        Array<long> mPopulation;
        Matrix<long> mMoleculeCount;
        
        // run statistics
        RunStatistics mRunStatistics;
    };
    
    	
    	
	inline void Simulation::SetNumSpecies(short num)
	{
		mSpeciesName.SetSize(num);
		mPopulation.SetSize(num);
		return;
	}
	
	
	
	inline short Simulation::NumSpecies() const
	{
		return mSpeciesName.Size();
	}
	
	
	
	inline void Simulation::SetSpeciesName(short index, const std::string &name)
	{
		mSpeciesName[index] = name;
        Capitalize(mSpeciesName[index]);
        
		return;
	}
    
    
    
    inline std::string Simulation::SpeciesName(short index) const
    {
        return mSpeciesName[index];
    }
    
    
    
    inline void Simulation::SetNumMaterials(short number)
    {
    	mMaterial.SetSize(number);
    	return;
    }
    
    
    
    inline Reaction* Simulation::AddReaction(const Reaction &reaction)
    {
    	mReactionList.push_back(reaction);
    	std::list<Reaction>::iterator iRL = --(mReactionList.end());
    	(*iRL).SetIterator(&*iRL);

    	return &*iRL;
    }
    


	inline void Simulation::SetReactionRate(Reaction *pReaction) const
	{
		short index = pReaction->ReactionIndex();
		short neighborNumber = pReaction->NeighborIndex();
		Element *pE = pReaction->GetElement();
	    	
		if (neighborNumber < 0) 
    		pReaction->SetRate(mChemicalSystem.Rate(index, pE));
		else 
    		pReaction->SetRate(DiffusionRate(pE, index, neighborNumber));

		return;
	}



	inline double Simulation::DiffusionRate(const Element* pE, short speciesNumber, short neighborNumber) const
	{
		short n = pE->GetOccupancyData()->NumOccupants(speciesNumber);
		if (n == 0)
			return 0.0;
	
   		double factor = pE->GetDiffusionFactor()->Factor(neighborNumber);
   		
    	double bigDe = mMaterial[pE->MaterialNumber()].DiffusionCoefficient(speciesNumber); 
    	
		double bigDn = mMaterial[pE->PN(neighborNumber)->MaterialNumber()].DiffusionCoefficient(speciesNumber);
		
		double bigD = 0.5 * (bigDe + bigDn);
		//double bigD = 2.0 * bigDe * bigDn / (bigDe + bigDn);
		
		return bigD * factor * n;
	}
}

#endif // _simulation_h_	
