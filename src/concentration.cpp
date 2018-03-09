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

#include <math.h>
#include <list>

#include "concentration.h"

using namespace NAMESPACE;
using namespace std;

double Concentration::Value(const Point *pPt, short materialNumber) const
{
	switch (mType) {
	case NO_FUNCTION_TYPE:
		ThrowException("Concentration::Value : unspecified function type");
		break;
			
	case UNIFORM:
		return Uniform(pPt, materialNumber);
		break;
			
	case GAUSSIAN:
		return Gaussian(pPt);
		break;
        
    case STEP_FUNCTION:
        return StepFunction(pPt);
        break;
		   
    case CONVERGENCE_TEST:
        return ConvergenceTest(pPt);
        break;  
       	   
    case DISCONTINUITY_TEST:
        return DiscontinuityTest(pPt);
        break;  
        
    case CUBIC_WAVE_EXACT:
        return CubicWaveExact(pPt);
        break;  
        
	default:
		ThrowException("Concentration::Value : unspecified function type");
		break;
	}
	
	return -1.0;
}



double Concentration::Uniform(const Point *pPt, short materialNumber) const
{
	// mParameter[0] is the concentration
	if (mParameter[0] < 0.0)
		ThrowException("Concentration::Uniform : concentration is negative");
	
	if ((short) mParameter[1] == materialNumber)
		return mParameter[0];
	else
		return 0.0;
}



double Concentration::Gaussian(const Point *pPt) const
{
	// Gaussian is of the form
	// 
	// A exp(-((x-x_0) / a)^2) * exp(-((y-y_0) / b)^2) * exp(-((z-z_0) / c)^2)
	//
	// mParameter[0] :	A
	// mParameter[1] : x_0
	// mParameter[2] : y_0
	// mParameter[3] : z_0
	// mParameter[4] : a
	// mParameter[5] : b
	// mParameter[6] : c
	
	if (mParameter[0] <= 0.0)
		ThrowException("Concentration::Gaussian : concentration is non-positive");
		
	if (mParameter[4] == 0.0)
		ThrowException("Concentration::Gaussian : zero X width");
		
	double xArg = (pPt->X() - mParameter[1]) / mParameter[4];
	double xTerm = exp(-xArg * xArg);

	if (mParameter[5] == 0.0)
		ThrowException("Concentration::Gaussian : zero Y width");
		
	double yArg = (pPt->Y() - mParameter[2]) / mParameter[5];
	double yTerm = exp(-yArg * yArg);
	
	if (mParameter[6] == 0.0)
		ThrowException("Concentration::Gaussian : zero Z width");
		
	double zArg = (pPt->Z() - mParameter[3]) / mParameter[6];
	double zTerm = exp(-zArg * zArg);
	
	return mParameter[0] * xTerm * yTerm * zTerm;
}



double Concentration::StepFunction(const Point *pPt) const
{
	// mParameter[0] is the coordinate direction
    // a positive sign means that the step function is non-zero to the right
    // of the step location, whereas a negative sign means it is non-zero
    // to the left of the step location
    long dir = NearestInteger(mParameter[0]);
    
    // mParameter[1] is the location of the step
    double locationOfPoint;
    switch (abs(dir)) {
    case 1:
        locationOfPoint = pPt->X();
        break;
        
    case 2:
        locationOfPoint = pPt->Y();
        break;
        
    case 3:
        locationOfPoint = pPt->Z();
        break;
    
    default:
        ThrowException("Concentration::StepFunction : bad coordinate direction");
        break;
    }
    
    if (dir < 0.0) {
    	if (locationOfPoint < mParameter[1])
        	return mParameter[2];
    	else
        	return 0.0;
    }
    else {
    	if (locationOfPoint > mParameter[1])
        	return mParameter[2];
    	else
        	return 0.0;
    }
}



double Concentration::ConvergenceTest(const Point *pPt) const
{
	long n = NearestInteger(mParameter[1]);
	long m = NearestInteger(mParameter[2]);
	long p = NearestInteger(mParameter[3]);
	
    double x = pPt->X();
    double y = pPt->Y();
    double z = pPt->Z();
    
    double term = -PI * PI * (n * n + m * m + p * p);
    term = exp(term * mParameter[4]);

	return mParameter[0] * (1.0 + 0.5 * term * cos(n * PI * x) * cos(m * PI * y) * cos(p * PI * z));
}



double Concentration::DiscontinuityTest(const Point *pPt) const
{
	static const double lambda = 9.0064951747691;
	static const double rho = 2.0651071015366;
	static const double sqrtFive = sqrt(5.0);
	
    double x = pPt->X();
    double bigX;
	if (x < 0.5) 
		bigX = rho * cos(lambda * x);
	else 
		bigX = cos(lambda * (1.0 - x) / sqrtFive);
	
    double term = exp(-lambda * lambda * mParameter[1]);

	return mParameter[0] * (1.0 + 0.25 * term * bigX);
}



double Concentration::CubicWaveExact(const Point *pPt) const
{
	double a0 = mParameter[0];
	double z0 = mParameter[1];
	double tau = mParameter[3];
	
	double z = pPt->X() - (tau / SQRT_TWO);
	double f = exp(-(z - z0) / SQRT_TWO);
	double beta = f / (1.0 + f);
	
	if (NearestInteger(mParameter[2]) == 2) {
		return a0 * beta;
	}
	else {
		return a0 - a0 * beta;
	}
}



void Concentration::SetType(const string &name)
{
	SetType(GetConcentrationType(name));
    return;
}



ConcentrationFunctionType Concentration::GetConcentrationType(const string &name)
{
	string tmp = name;
	LowerCase(tmp);
	
	if (tmp == "uniform")
		return UNIFORM;
    
	if (tmp == "gaussian")
		return GAUSSIAN;
    
    if (tmp == "deltafunction")
        return DELTA_FUNCTION;
    
    if (tmp == "stepfunction")
        return STEP_FUNCTION;
    
    if (tmp == "convergencetest")
        return CONVERGENCE_TEST;
     
    if (tmp == "discontinuitytest")
        return DISCONTINUITY_TEST;
    
    if (tmp == "cubicwaveexact")
        return CUBIC_WAVE_EXACT;
    
	ThrowException("GetConcentrationType : bad input string : " + name);
    
	return NO_FUNCTION_TYPE;
}
