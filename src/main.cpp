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
#include "parser.h"
#include "realmatrix.h"
#include "meshgenerator.h"
#include "clock.h"

#include <iostream>
#include <fstream>
#include "utility.h"

//#include "Profiler.h"

using namespace NAMESPACE;
using namespace std; 

void MakeMeshFile(void);
void ConvertMeshFile(void);
void ReadPolygonFile(void);
double Ran013(void);

int main()
{	
	#if __profile__
		if (!ProfilerInit(collectDetailed, bestTimeBase, 2000, 500))
		{
	#endif
	
    try { 	
    	//RandomNumberTest();
    	
        //MakeMeshFile();
        
        //PruneTimeAverageFile();
        
        //Test();
        
        //ReadPolygonFile();
        
        Simulation sim;
        sim.Run("/Users/dave/Projects/ReDi/SimulationFiles/bonner1.txt");
    }
    catch (exception &standardException) {
        HandleException(standardException);
    }
    catch (string &message) {
        HandleException(message);
    }
	
			
	#if __profile__
		ProfilerDump("\pMain:Users:dave:Projects:ReDi:ReDi.prof");
		ProfilerTerm();
		}
	#endif

    return 0;
}



void MakeMeshFile()
{
    MeshGenerator mesh;
    
    //mesh.MeshSphere(1.0, 5);
    //mesh.MeshPill(0.5, 2.0, 12);
	//mesh.MeshRectangularSolidWithBricks(5.0, 5.0, 5.0, 20);
	mesh.MeshSphereWithBricks(6.0, 1.0);
	//mesh.MeshRectangle(1.0, 1.0, 10);
	//mesh.MeshLineEvenly(0.0, 1.0, 22);
	//mesh.MeshLineRandomly(0.0, 1.0, 40);

	//mesh.CheckAll();
	
	mesh.MaterialNumberTest();
    //display.ShowMesh(mesh);
    
    //display.SaveCurrentImageAsJPEG("/Users/dave/Projects/pill.jpg");
    
    //cout << mesh.NumElements() << endl;
    
    //mesh.WriteMatlabEdgeFile();
    mesh.WriteToFile("/Users/dave/Projects/ReDi/MeshFiles/sphere_5_1.dat");
    return;
}



void ConvertMeshFile()
{
	ifstream input("/Users/dave/Projects/ReDi/MeshFiles/wavetube.dat");
	ofstream output("/Users/dave/Projects/ReDi/MeshFiles/newmesh.dat");
	
	long numVertices;
	input >> numVertices;
	output << numVertices << endl;
	double x, y, z;
	long id;
	for (long i = 0; i < numVertices; ++i) {
		input >> id;
		input >> x;
		input >> y;
		input >> z;
		
		output << id << " ";
		output << x << " ";
		output << y << " ";
		output << z << endl;
	}
	
	// read elements
	long numElements;
	input >> numElements;
	output << numElements << endl;
	
	Array<short> eType(numElements);
	
	long eID;
	for (long i = 0; i < numElements; ++i) {
		input >> eID;
		input >> eType[i];
		
		output << eID << " ";
		output << eType[i] << " ";
		output << 0 << " ";
		
		short numVertices = (eType[i] == TRIANGLE) ? 3 : 4;
		if (eType[i] == BLOB)
			numVertices = 0;
		
		for (short j = 0; j < numVertices; ++j) {
			input >> id;
			output << id  << " ";
		}
		output << endl;
	}

	// read element neighbors
	for (long i = 0; i < numElements; ++i) {
		short numNeighbors = (eType[i] == TRIANGLE) ? 3 : 4;
		if (eType[i] == BLOB)
			numNeighbors = 0;
		
		for (short j = 0; j < numNeighbors; ++j) {
			input >> id;
			output << id << " ";
		}
		output << endl;
	}
	

	return;
}



void ReadPolygonFile()
{
	MeshGenerator meshGen;
	meshGen.MakeMeshFromPolygonFile("/Users/dave/Projects/ReDi/v400.dat");
	
	meshGen.WriteToFile("/Users/dave/Projects/ReDi/MeshFiles/voronoi_400.dat");

	return;
}




unsigned long x = (1 ^  521288629) & 0xffffffff;
unsigned long y = (1 ^  362436069) & 0xffffffff;
unsigned long z = (1 ^   16163801) & 0xffffffff;
unsigned long n = (1 ^ 1131199209) & 0xffffffff;
unsigned long c = y > z;

double Ran013()
{
  // The mask 0xffffffff is neccessary in some places to assure that arithmetics
  // is performed modulo 2^32 to make the generater portable to any word length 
  // larger than 2^32.

  long s;
  if (y > x+c) { 
    s = y - (x+c); 
    c = 0; 
  } else {
    s = (y - (x+c) - 18) & 0xffffffff; // mask is neccessary here
    c = 1;
  }
  x = y;
  y = z;   

	//return (((z = s) + (n = 69069 * n + 1013904243)) & 0xffffffff);
  return (((z = s) + (n = 69069 * n + 1013904243)) & 0xffffffff)*2.328306437e-10;  // and here
}
