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

#ifndef _redienums_h_
#define _redienums_h_

#include "namespace.h"

namespace NAMESPACE {
	enum StopFlag{NO_STOP, 
				  END_TIME_REACHED, 
				  NO_REACTION_POSSIBLE, 
				  CONTACT_WITH_BOUNDARY, 
				  ZERO_POPULATION};
				  
	enum ConcentrationFunctionType{NO_FUNCTION_TYPE, 
                                   UNIFORM, 
                                   GAUSSIAN, 
                                   DELTA_FUNCTION, 
                                   STEP_FUNCTION,
                                   CONVERGENCE_TEST,
                                   DISCONTINUITY_TEST,
                                   CUBIC_WAVE_EXACT};
    
    // note that in this enum ENSEMBLE_DATA_OUTPUT *must* be last
    enum DataOutputType{CONSOLE_DATA_OUTPUT, 
    					MESH_DATA_OUTPUT, 
    					TIME_DATA_OUTPUT, 
    					XYDATA_OUTPUT, 
    					TIME_AVERAGE,
    					TIME_AVERAGE_OUTPUT,
    					ENSEMBLE_AVERAGE_TIME_OUTPUT,
    					ENSEMBLE_DATA_OUTPUT};
}

#endif // _redienums_h_
