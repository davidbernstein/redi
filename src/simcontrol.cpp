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

#include "simcontrol.h"

using namespace std;
using namespace NAMESPACE;

bool SimControl::NeedRefinementData()
{
	return (mInitialRefinementLevel > 0) || mMeshRefinementOn || mRefinementTestOn || mUnRefinementTestOn;
}



short SimControl::ComputeNumDataOutputTypes() const
{
	short i = 0;
	DataOutputType type = (DataOutputType) i;
	while (type != ENSEMBLE_DATA_OUTPUT) {
		type = (DataOutputType) i;
		++i;
	}
		
	return i;
}
