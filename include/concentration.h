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

#ifndef _concentration_h_
#define _concentration_h_

#include "redienums.h"
#include "rediconst.h"

#include "point.h"

namespace NAMESPACE {
	class Concentration {
	public:
		// Constructor
		Concentration(ConcentrationFunctionType type = NO_FUNCTION_TYPE);

		// Destructor
		~Concentration(void) { };

		// value
		double Value(const Point *pPt, short materialNumber = 0) const;
		
		// type
		void SetType(ConcentrationFunctionType type);
		void SetType(const std::string &name);
        ConcentrationFunctionType Type(void) const;
        ConcentrationFunctionType GetConcentrationType(const std::string &name);
		
		// number of parameters
		short NumParameters(void);
		void SetParameter(short index, double p);
        double Parameter(short index) const;
		
	private:
		double Gaussian(const Point *pPt) const;
		double Uniform(const Point *pPt, short materialNumber = 0) const;
        double StepFunction(const Point *pPt) const;
        double ConvergenceTest(const Point *pPt) const;
        double DiscontinuityTest(const Point *pPt) const;
		double CubicWaveExact(const Point *pPt) const;
		
	private:
		// type of function
		ConcentrationFunctionType mType;
		
		// parameters
		Array<double> mParameter;
	};
	
	
	
	inline Concentration::Concentration(ConcentrationFunctionType type)
	{
		SetType(type);
		return;
	}
	
	
	
	inline void Concentration::SetType(ConcentrationFunctionType type)
	{
		if (type == NO_FUNCTION_TYPE) {
			mParameter.SetSize(0);
			mType = NO_FUNCTION_TYPE;
			return;
		}
		
		if (type == UNIFORM) {
			mParameter.SetSize(2);
			mType = UNIFORM;
			return;
		}
			
		if (type == GAUSSIAN) {
			mParameter.SetSize(7);
			mType = GAUSSIAN;
			return;
		}
        
		if (type == DELTA_FUNCTION) {
			mParameter.SetSize(4);
			mType = DELTA_FUNCTION;
			return;
		}
        
		if (type == STEP_FUNCTION) {
			mParameter.SetSize(3);
			mType = STEP_FUNCTION;
			return;
		}
		
		if (type == CONVERGENCE_TEST) {
			mParameter.SetSize(5);
			mType = CONVERGENCE_TEST;
			return;
		}
        
		if (type == DISCONTINUITY_TEST) {
			mParameter.SetSize(2);
			mType = DISCONTINUITY_TEST;
			return;
		}
		
		if (type == CUBIC_WAVE_EXACT) {
			mParameter.SetSize(4);
			mType = CUBIC_WAVE_EXACT;
			return;
		}
		
		ThrowException("Concentration::SetType : bad type");
		
		return;
	}
	
	
	
	inline short Concentration::NumParameters()
	{
		return mParameter.Size();
	}
	
	
	
	inline void Concentration::SetParameter(short index, double p)
	{
		mParameter[index] = p;
		return;
	}
    
    
    
    inline double Concentration::Parameter(short index) const
    {
        return mParameter[index];
    }
    
    
    
    inline ConcentrationFunctionType Concentration::Type() const
    {
        return mType;
    }
}


#endif // _concentration_h_

