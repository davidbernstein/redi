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

#include "chemicalsystem.h"

#include <iostream>

using namespace NAMESPACE;
using namespace std;

double ChemicalSystem::Rate(short reactionIndex, const Element *pE) const
{
	Occupancy *pOcc = pE->GetOccupancyData();
	
	double rate = RateConstant(reactionIndex);
	for (short i = 0; i < mNumSpecies; ++i) {
		short coef = ReactantCoefficient(reactionIndex, i);
		short n = pOcc->NumOccupants(i);
		
		if (n < coef)
			return 0.0;
		
		rate *= (coef > 0) ? BigCombination(n, coef) : 1;
	}
	
	short n = SumReactantCoefficients(reactionIndex) - 1;
	
	rate /= Power(pE->Volume(), n);
	
	//return rate * mRateFactor[reactionIndex];
	rate *= mRateFactor[reactionIndex];
	return rate;
}



bool ChemicalSystem::ReactionDependentOnReaction(short indexA, short indexB) const
{
	for (short i = 0; i < NumSpecies(); ++i) {
		if (ReactantCoefficient(indexB, i) > 0) {
			if (ReactionDependentOnSpecies(indexA, i))
				return true;
		}
		
		if (ProductCoefficient(indexB, i) > 0) {
			if (ReactionDependentOnSpecies(indexA, i))
				return true;
		} 
	}
	
	return false;
}



bool ChemicalSystem::ReactionDependentOnSpecies(short reactionIndex, short speciesNumber) const
{
	return ReactantCoefficient(reactionIndex, speciesNumber) > 0;
}



void ChemicalSystem::ComputeRateFactors()
{
	mRateFactor.SetSize(mNumReactions);
	
	for (short r = 0; r < mNumReactions; ++r) {
		double factor = 1;
		for (short i = 0; i < mNumSpecies; ++i) 
			factor *= SmallFactorial(mReactantCoefficient(r, i));
	
		mRateFactor[r] = factor;
	}
	
	return;
}



void ChemicalSystem::SetSpeciesNames(const Array<std::string> &names)
{
	mSpeciesName = names;
	return;
}



short ChemicalSystem::WhichSpecies(std::string speciesName)
{
	for (short i = 0; i < mSpeciesName.Size(); ++i) {
		if (speciesName == mSpeciesName[i])
			return i;
	}
	
	ThrowException("ChemicalSystem::WhichSpecies : name not found: " + speciesName);
	
	return -1;
}



void ChemicalSystem::SetReactionMaterialActivity(short index, const Array<bool> &activity)
{
	if (activity.Size() != mNumReactions)
		ThrowException("ChemicalSystem::SetReactionMaterialActivity : activity array wrong size");
	
	for (short i = 0; i < mNumReactions; ++i)
		mReactionMaterialActivity(index, i) = activity[i];
	
	return;
}



void ChemicalSystem::InitializeReactionMaterialActivity(short numMaterials)
{
	mReactionMaterialActivity.SetSize(mNumSpecies, numMaterials);
	for (short i = 0; i < mNumSpecies; ++i) {
		for (short j = 0; j < numMaterials; ++j)
			mReactionMaterialActivity(i,j) = true;
	}
	
	return;
}



void ChemicalSystem::PrintReactions() const
{
	if (mNumReactions == 0)
		return;
	
	cout << "chemical system:" << endl;
	for (short r = 0; r < mNumReactions; ++r) {
		PrintReactionLine(r, mReactantCoefficient);
		cout << " --> ";
		PrintReactionLine(r, mProductCoefficient);
		cout << endl;
	}
	

	return;
}



void ChemicalSystem::PrintReactionLine(short reactionIndex, const Matrix<short> &coef) const
{
	short pos = LastNonZeroPosition(reactionIndex, coef);
	
	short count = 0;
	for (short i = 0; i <= pos; ++i) {
		if (coef(reactionIndex, i) > 0) {
			if (coef(reactionIndex, i) == 1) 
				cout << mSpeciesName[i];
			
			if (coef(reactionIndex, i) > 1) 
				cout << coef(reactionIndex, i) << mSpeciesName[i];
			
			if (i < pos)
				cout << " + ";
			
			++count;
		}
	}
	
	if (count == 0)
		cout << "null";
	
	
	return;
}



short ChemicalSystem::LastNonZeroPosition(short reactionIndex, const Matrix<short> &coef) const
{
	short pos = 0;
	
	for (short i = 0; i < mNumSpecies; ++i) {
		if (coef(reactionIndex, i) > 0) {
			pos = i;
		}
	}
	
	return pos;
}