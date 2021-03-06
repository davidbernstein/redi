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

#include <fstream>
#include "mesh.h"

using namespace std;
using namespace NAMESPACE;


void Mesh::WriteMCFile(const string &fileName, Array<short> &elementGroupNumber, Array<Array<short> > &faceType) const
{
	ofstream file(fileName.c_str());
	file.setf(ios::scientific);

	WriteMCFile(file, elementGroupNumber, faceType);
	
	file.close();

	return;
}



void Mesh::WriteMCFile(ofstream &file, Array<short> &elementGroupNumber, Array<Array<short> > &faceType) const
{
	// write file header
	WriteMCHeader(file);
	
	// write vertices
	WriteMCVertices(file);
	
	// write elements
	WriteMCElements(file, elementGroupNumber, faceType);
	
	// write footer
	WriteMCFooter(file);
	
	return;
}



void Mesh::WriteMCHeader(ofstream &file) const
{
	file << "%" << endl;
	file << "mcsf_begin=1;" << endl;
	file << "dim=3;" << endl;
	file << "dimii=3;" << endl;
	file << "vertices=" << NumVertices() << ";" << endl;
	file << "simplices=" << mTetrahedronList.size() << ";" << endl;

	return;
}



void Mesh::WriteMCFooter(ofstream &file) const
{
	file << "mcsf_end=1;" << endl;
	return;
}



void Mesh::WriteMCVertices(ofstream &file) const
{
	file << "vert=[" << endl;
	list<Vertex>::const_iterator iVl;
	for (iVl = mVertexList.begin(); iVl != mVertexList.end(); ++iVl) {
		file << (*iVl).ID() << " " << "0 ";
		file << (*iVl).X() << " " << (*iVl).Y() << " " << (*iVl).Z();
		file << endl;
	}
	file << "];" << endl;

	return;
}



void Mesh::WriteMCElements(ofstream &file, Array<short> &elementGroupNumber, Array<Array<short> > &faceType) const
{
	file << "simp=[" << endl;
	
	list<Tetrahedron>::const_iterator iTe;
	for (iTe = mTetrahedronList.begin(); iTe != mTetrahedronList.end(); ++iTe) {
		// write simplex id
		long id = (*iTe).ID();
		file << id << " ";

		// write simplex group number
		file << elementGroupNumber[id] << " ";

		// material number is zero for all simplices
		file << "0 ";
		
		// write simplex face type
		for (short i = 0; i < 4; ++i)
			file << faceType[id][i] << " ";
		
		// write vertices
		(*iTe).WriteVertexIDsToFile(file);
		
		file << endl;
	}
	
	file << "];" << endl;
	
	return;
}
