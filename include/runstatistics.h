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

#ifndef _runstatistics_h_
#define _runstatistics_h_

#include "redienums.h"
#include "clock.h"

namespace NAMESPACE {
	class RunStatistics {
    public:
		RunStatistics(void);
        ~RunStatistics(void) { };
		
		void Initialize(void);
		
		// run time
		double AverageRunTime(void) const;
		
		// number of runs
		void IncrementNumberOfRuns(void);
		
		// number of steps
		void AddToTotalNumberOfSteps(long numSteps);
		double AverageNumberOfSteps(void) const;
		
		// clock
		void StartClock(void);
		void StopClockAndPrintTime(void);
		
    private:
    	// total run time 
    	double mRunTime;
    	
    	// total number of steps
    	double mNumSteps;
    	
    	// number of runs
    	long mNumRuns;
    	
    	// clock for timing
        Clock mClock;
	};
    
    
    
    inline RunStatistics::RunStatistics() 
    {
    	Initialize();
        return;
    }
    
    
    
    inline void RunStatistics::Initialize()
    {
    	mRunTime = 0.0;
        mNumSteps = 0;
        mNumRuns = 0;
    	
    	mClock.Stop();
    	
        return;
    }
    
    
    inline double RunStatistics::AverageRunTime() const
    {
    	if (mNumRuns == 0)
    		ThrowException("RunStatistics::AverageRunTime : number of runs is zero");
    	
    	return mRunTime / mNumSteps;
    }
		
	
	
	inline void RunStatistics::IncrementNumberOfRuns()
	{
		++mNumRuns;
		return;
	}
		
		
	
	inline void RunStatistics::AddToTotalNumberOfSteps(long numSteps)
	{
		mNumSteps += numSteps;
		return;
	}
	
	
	
	inline double RunStatistics::AverageNumberOfSteps() const
	{
		if (mNumRuns == 0)
    		ThrowException("RunStatistics::AverageRunTime : number of runs is zero");
    	
    	return mNumSteps / mNumSteps;
    }
    
    
		
	inline void RunStatistics::StartClock()
	{
		mClock.Start();
		return;
	}
	
	
	
	inline void RunStatistics::StopClockAndPrintTime()
	{
		mRunTime += mClock.StopAndPrintTime("Run time ");
		return;
	}
}


#endif // _runstatistics_h_


