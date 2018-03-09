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

#ifndef _chemicalsystem_h_
#define _chemicalsystem_h_

#include <string>
#include "array.h"
#include "matrix.h"
#include "element.h"

namespace NAMESPACE {
	class ChemicalSystem
	{
	public:
		// Constructor
		ChemicalSystem(void) { };

		// Destructor
		~ChemicalSystem(void) { };

		// number of species
		void SetNumReactionsAndSpecies(short numReactions, short numSpecies);
		short NumSpecies(void) const;
		short NumReactions(void) const;
		
		// rate
		void SetRateConstant(short reactionIndex, double rate);
		double RateConstant(short reactionIndex) const;
		double Rate(short reactionIndex, const Element *pE) const;
        
		// coefficients
		void SetReactantCoefficients(short reactionIndex, const Array<short> &coef);
		void SetProductCoefficients(short reactionIndex, const Array<short> &coef);
		short ReactantCoefficient(short reactionIndex, short speciesIndex) const;
		short ProductCoefficient(short reactionIndex, short speciesIndex) const;
		short SumReactantCoefficients(short reactionIndex) const;
		
		// dependent reactions
		bool ReactionDependentOnReaction(short indexA, short indexB) const;
		bool ReactionDependentOnSpecies(short reactionIndex, short speciesNumber) const;
		
		// material activity
		void InitializeReactionMaterialActivity(short numMaterials);
		bool ReactionActive(short index, short materialNumber) const;
		void SetReactionMaterialActivity(short index, const Array<bool> &activity);
		
		// rate factors
		void ComputeRateFactors(void);
		double RateFactor(short reactionIndex) const;
		
		// names
		void SetSpeciesNames(const Array<std::string> &names);
		short WhichSpecies(std::string speciesName);

		// IO
		void PrintReactions(void) const;
		
	private:
		void PrintReactionLine(short reactionIndex, const Matrix<short> &coef) const;
		short LastNonZeroPosition(short reactionIndex, const Matrix<short> &coef) const;

	
	private:
		// number of species
		short mNumSpecies;
		
		// number of reactions
		short mNumReactions;
		
		// rate constants
		Array<double> mRateConstant;
		
		// reactant coefficients
		// mReactantCoefficient(i, j): i = which reaction, j = which species
		Matrix<short> mReactantCoefficient;
		
		// product coefficients
		Matrix<short> mProductCoefficient;
	 
	 	// which reactions are active in which materials
	 	Matrix<bool> mReactionMaterialActivity;
	 	
        // species names
        Array<std::string> mSpeciesName;
        
        // factor used in rate computation
        Array<double> mRateFactor;
	};



	inline short ChemicalSystem::NumSpecies() const
	{
		return mNumSpecies;
	}
	
	
	
	inline short ChemicalSystem::NumReactions() const
	{
		return mNumReactions;
	}
	

	
	inline void ChemicalSystem::SetNumReactionsAndSpecies(short numReactions, short numSpecies)
	{
		mNumSpecies = numSpecies;
		mNumReactions = numReactions;
		
		mRateConstant.SetSize(numReactions);
		
		mReactantCoefficient.SetSize(numReactions, numSpecies);
		mProductCoefficient.SetSize(numReactions, numSpecies);
		
		for (short i = 0; i < numReactions; ++i) {
			for (short j = 0; j < numSpecies; ++j) {
				mReactantCoefficient(i, j) = 0;
				mProductCoefficient(i, j) = 0;
			}
		}
		
		
		return;
	}
	
	
	
	inline void ChemicalSystem::SetReactantCoefficients(short reactionIndex, const Array<short> &coef)
	{
        for (short i = 0; i < coef.Size(); ++i)
            mReactantCoefficient(reactionIndex, i) = coef[i];
            
		return;
	}
	
	
	
	inline void ChemicalSystem::SetProductCoefficients(short reactionIndex, const Array<short> &coef)
	{
        for (short i = 0; i < coef.Size(); ++i)
            mProductCoefficient(reactionIndex, i) = coef[i];

		return;
	}
	
	
	
	inline short ChemicalSystem::ReactantCoefficient(short reactionIndex, short speciesIndex) const
	{
		return mReactantCoefficient(reactionIndex, speciesIndex);
	}
	
	
	
	inline short ChemicalSystem::SumReactantCoefficients(short reactionIndex) const
	{
		short sum = 0;
		for (short i = 0; i < mNumSpecies; ++i) 
			sum += mReactantCoefficient(reactionIndex, i);
		
		return sum;
	}
	
	
	
	inline short ChemicalSystem::ProductCoefficient(short reactionIndex, short speciesIndex) const
	{
		return mProductCoefficient(reactionIndex, speciesIndex);
	}
	
	
	
	inline void ChemicalSystem::SetRateConstant(short reactionIndex, double rate)
	{
		mRateConstant[reactionIndex] = rate;
		return;
	}
	
	
	
	inline double ChemicalSystem::RateConstant(short reactionIndex) const
	{
		return mRateConstant[reactionIndex];
	}
	
	
	
	inline double ChemicalSystem::RateFactor(short reactionIndex) const
	{
		return mRateFactor[reactionIndex];
	}
	
	
	
	inline bool ChemicalSystem::ReactionActive(short index, short materialNumber) const
	{
		return mReactionMaterialActivity(index, materialNumber);
	}
}


#endif // _chemicalsystem_h_	
