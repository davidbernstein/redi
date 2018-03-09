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

#ifndef _simcontrol_h_
#define _simcontrol_h_

#include "constants.h"
#include "redienums.h"

namespace NAMESPACE {
	class SimControl {
    public:
		SimControl(void);
        ~SimControl(void) { };
		
        void SetRunTime(double startTime, double endTime);
        double EndTime(void) const;
        double StartTime(void) const;
        
        void SetNumberOfRuns(long numRuns);
        long NumberOfRuns(void) const;
        
        void SetToStopOnContactWithBoundary(void);
        bool StopOnContactWithBoundary(void) const;
        
        void SetStopOnZeroPopulationSize(short numSpecies);
        void SetToStopOnZeroPopulation(short speciesIndex);
        bool ZeroPopulationStop(void) const;
        bool ZeroPopulationStop(short speciesIndex) const;
        
        void SetOutputInterval(double value, DataOutputType dataType);
        double OutputInterval(DataOutputType dataType) const;
        bool Output(double time, double &outputTime, DataOutputType dataType);
        bool OutputOn(DataOutputType dataType);
        void ResetNextOutputTimes(void);
        short NumDataOutputTypes(void) const;
        bool OutputEnsembleData(long n);
        
        short InitialRefinementLevel(void) const;
        void SetInitialRefinementLevel(short value);
        
        bool MeshRefinementOn(void) const;
        void TurnOnMeshRefinement(void);
        bool RefinementTestOn(void) const;
        void TurnOnRefinementTest(void);
        bool UnRefinementTestOn(void) const;
        void TurnOnUnRefinementTest(void);
        bool NeedRefinementData(void);
        
        void TurnOnConvergenceTestError(void);
        bool ConvergenceTestErrorOn(void) const;
        
        void TurnOnRandomWalkTest(void);
        bool RandomWalkTestOn(void) const;
        
        void TurnOnCubicTestError(void);
        bool CubicTestErrorOn(void) const;
        
        void TurnOnEndTimeOutput(void);
        bool OutputEndTime(void) const;
        
        short SpatialAverageNeighborhoodSize(void) const;
        void SetSpatialAverageNeighborhoodSize(short value);
        
    private:
    	short ComputeNumDataOutputTypes(void) const;
    	
    private:
        // time of run
        double mStartTime;
		double mEndTime;
        
        // number of runs in ensemble
        long mNumberOfRuns;
        
        // stop run if any molecule contacts a boundary
        bool mStopOnContactWithBoundary;
        
        // stop on population count
        bool mZeroPopulationStop;
        Array<bool> mStopOnZeroPopulation;
        
        // mesh refinement on 
        bool mMeshRefinementOn;
        
        // maximum refinement level
        long mMaxRefinementLevel;
        
        // initial refinement level
        short mInitialRefinementLevel;
        
        // refinement testing
        bool mRefinementTestOn;
        bool mUnRefinementTestOn;
        
        // interval between writing data
        short mNumDataOutputTypes;
        Array<double> mOutputInterval;
        Array<double> mNextOutputTime;
        Array<bool> mPerformOutput;
        
        // error output
        bool mConvergenceTestErrorOn;
        bool mRandomWalkTestOn;
        bool mCubicTestErrorOn;
        bool mOutputEndTime;
        
        // time averaging
        long mNumberOfTimeAveragesPerOutput;
        
        // spatial averaging
        short mSpatialAverageNeighborhoodSize;
	};
    
    
    
    inline SimControl::SimControl() 
    {
        mStartTime = 0.0;
        mEndTime = 0.0;
        
        mNumberOfRuns = 1;
                
        mStopOnContactWithBoundary = false;
        mMaxRefinementLevel = MACHINE_LARGEST_LONG;
        
        mInitialRefinementLevel = 0;
        mMeshRefinementOn = false;
        mRefinementTestOn = false;
        mUnRefinementTestOn = false;
        
        mNumDataOutputTypes = ComputeNumDataOutputTypes();
        mOutputInterval.SetSize(mNumDataOutputTypes);
        mNextOutputTime.SetSize(mNumDataOutputTypes);
        mPerformOutput.SetSize(mNumDataOutputTypes);
        
        for (short i = 0; i < mNumDataOutputTypes; ++i) {
        	mOutputInterval[i] = MACHINE_LARGEST_DOUBLE;
        	mNextOutputTime[i] = mStartTime;
        	mPerformOutput[i] = false;
        }
        
        mConvergenceTestErrorOn = false;
        mRandomWalkTestOn = false;
        mCubicTestErrorOn = false;
        mOutputEndTime = false;
        
        mNumberOfTimeAveragesPerOutput = 0;
        
        mSpatialAverageNeighborhoodSize = 0;
        
        mZeroPopulationStop = false;
        
            	
        return;
    }



    inline void SimControl::SetRunTime(double startTime, double endTime)
    {
        if (startTime > endTime)
            ThrowException("Simulation::SetRunTime : start time greater than end time");
    
        mStartTime = startTime;
        mEndTime = endTime;
				
        return;
    }
    
    
    
    inline double SimControl::StartTime() const
    {
        return mStartTime;
    }
    
    
    
    inline double SimControl::EndTime() const
    {
        return mEndTime;
    }
    
    
    
    inline void SimControl::SetNumberOfRuns(long numRuns)
    {
        if (numRuns < 1)
            ThrowException("SimControl::SetNumberOfRuns : number not positive");
            
        mNumberOfRuns = numRuns;
        return;
    }
    
    
    
    inline long SimControl::NumberOfRuns() const
    {
        return mNumberOfRuns;
    }
        
    
    
    inline void SimControl::SetToStopOnContactWithBoundary()
    {
        mStopOnContactWithBoundary = true;
        return;
    }
    
    
    
    inline bool SimControl::StopOnContactWithBoundary() const
    {
        return mStopOnContactWithBoundary;
    }
    
    
    
    inline short SimControl::NumDataOutputTypes() const
    {
    	return mNumDataOutputTypes;
    }
    
    
    
    inline void SimControl::SetOutputInterval(double value, DataOutputType dataType)
    {
        if (value <= 0.0)
            ThrowException("SimControl::SetOutputInterval : non-positive value");
        
        short index = (short) dataType;
    
        mOutputInterval[index] = value;
        mNextOutputTime[index] = mStartTime + value;
        mPerformOutput[index] = true;
        
        return;
    }
    
    
    
    inline double SimControl::OutputInterval(DataOutputType dataType) const
    {
    	return mOutputInterval[(short) dataType];
    }
    
    
    
    inline void SimControl::ResetNextOutputTimes()
    {
    	for (short i = 0; i < mNumDataOutputTypes; ++i) {
    		if (i != (short) ENSEMBLE_DATA_OUTPUT) {
    			mNextOutputTime[i] = mStartTime;
        		//mNextOutputTime[i] = mStartTime + mOutputInterval[i];
    		}
        }
        
        return;
    }
    
    
    
    inline bool SimControl::Output(double time, double &outputTime, DataOutputType dataType)
    {
    	short index = (short) dataType;
    	
    	if (mPerformOutput[index] == false)
    		return false;
    	
        if (time < mNextOutputTime[index]) {
            return false;
        }
        else {
        	outputTime = mNextOutputTime[index];
            mNextOutputTime[index] += mOutputInterval[index];
            return true;
        }
    }
    
    
    
    inline bool SimControl::OutputOn(DataOutputType dataType)
    {
    	return mPerformOutput[(short) dataType];
    }
    
    
    
    inline bool SimControl::OutputEnsembleData(long n)
    {
    	short index = (short) ENSEMBLE_DATA_OUTPUT;
    	
    	if (mPerformOutput[index] == false)
    		return false;
    	
    	if (mOutputInterval[index] >= 1.0) {
    		if (n < mNextOutputTime[index]) {
    			return false;
    		}
    		else {
    			mNextOutputTime[index] += mOutputInterval[index];
           		return true;
           	}
    	}
    	else {
    		double logn = log10((double) n);
    		if (logn < mNextOutputTime[index]) {
    			return false;
    		}
    		else {
    			mNextOutputTime[index] += mOutputInterval[index];
    			
    			while (logn > mNextOutputTime[index])
    				mNextOutputTime[index] += mOutputInterval[index];

           		return true;
           	}
    	}
    	
    	
    	return false;
    }
    
    
    
    inline short SimControl::InitialRefinementLevel() const
    {
    	return mInitialRefinementLevel;
    }
    
    
    
    inline void SimControl::SetInitialRefinementLevel(short value)
    {
    	if (value <= 0.0)
            return;
        
    	mInitialRefinementLevel = value;
    	
    	return;
    }
    
    
    
    inline bool SimControl::MeshRefinementOn() const
    {
    	return mMeshRefinementOn;
    }
    
    
    
    inline void SimControl::TurnOnMeshRefinement()
    {
    	mMeshRefinementOn = true;
    	return;
    }
    
    
    
    inline bool SimControl::RefinementTestOn() const
    {
    	return mRefinementTestOn;
    }

    
    
    inline void SimControl::TurnOnRefinementTest()
    {
    	mRefinementTestOn = true;
    	return;
    }
        
    
    
    inline bool SimControl::UnRefinementTestOn() const
    {
    	return mUnRefinementTestOn;
    }

    
    
    inline void SimControl::TurnOnUnRefinementTest()
    {
    	mUnRefinementTestOn = true;
    	return;
    }
    
    
    
    inline void SimControl::TurnOnConvergenceTestError()
    {
    	mConvergenceTestErrorOn = true;
    	return;
    }
    
    
    
    inline bool SimControl::ConvergenceTestErrorOn() const
    {
    	return mConvergenceTestErrorOn;
    }
    
    
    
    inline void SimControl::TurnOnCubicTestError()
    {
    	mCubicTestErrorOn = true;
    	return;
    }
    
    
    
    inline bool SimControl::CubicTestErrorOn() const
    {
    	return mCubicTestErrorOn;
    }
    
    
    
    inline void SimControl::TurnOnRandomWalkTest()
    {
    	mRandomWalkTestOn = true;
    	return;
    }
    
    
    
    inline bool SimControl::RandomWalkTestOn() const
    {
    	return mRandomWalkTestOn;
    }
    
    
    
    inline void SimControl::TurnOnEndTimeOutput()
    {
    	mOutputEndTime = true;
    	return;
    }
    
    
    
    inline bool SimControl::OutputEndTime() const
    {
    	return mOutputEndTime;
    }
    
    
      
    inline short SimControl::SpatialAverageNeighborhoodSize() const
    {
    	return mSpatialAverageNeighborhoodSize;
    }
    
    
    
    inline void SimControl::SetSpatialAverageNeighborhoodSize(short value)
    {	
    	mSpatialAverageNeighborhoodSize = (value < 0) ? 0 : value;
    }
    
    
    
    inline void SimControl::SetStopOnZeroPopulationSize(short numSpecies)
    {
    	mStopOnZeroPopulation.SetSize(numSpecies);
    	
    	for (short i = 0; i < numSpecies; ++i)
    		mStopOnZeroPopulation[i] = false;
    }
    
    
    
   	inline void SimControl::SetToStopOnZeroPopulation(short speciesIndex)
   	{
   		mStopOnZeroPopulation[speciesIndex] = true;
   		mZeroPopulationStop = true;
   		return;
   	}
   	
   	
   	
   	inline bool SimControl::ZeroPopulationStop(short speciesIndex) const
   	{
   		return mStopOnZeroPopulation[speciesIndex];
   	}
   	
   	
   	
   	inline bool SimControl::ZeroPopulationStop() const
   	{
   		return mZeroPopulationStop;
   	}
}


#endif // _simcontrol_h_


