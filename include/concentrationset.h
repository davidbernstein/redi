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

#ifndef _concentrationset_h_
#define _concentrationset_h_

#include "concentration.h"

namespace NAMESPACE {
	class Element;
}

namespace NAMESPACE {
	class ConcentrationSet {
	public:
		// Constructor
		ConcentrationSet() { };

		// Destructor
		~ConcentrationSet(void) { };

		// operators
		Concentration& operator[](long index);
		const Concentration& operator[](long index) const;
		
		// size
		void SetSize(short size);
		short Size(void) const;
		
		// delta function centroid
		void SetDeltaFunctionElement(short index, const Element *pElement);
		const Element* DeltaFunctionElement(short index) const;
		
	private:
		// concentrations
		Array<Concentration> mConcentration;
		
		// stuff having to do with delta functions
		Array<Element*> mDeltaFunctionElement;
	};
	
	
	
	inline void ConcentrationSet::SetSize(short size)
	{
		mConcentration.SetSize(size);
		mDeltaFunctionElement.SetSize(size);
		
		for (short i = 0; i < size; ++i)
			mDeltaFunctionElement[i] = NULL;
		
		return;
	}
	
	
	
	inline short ConcentrationSet::Size() const
	{
		return mConcentration.Size();
	}
	
	
	
	inline Concentration& ConcentrationSet::operator[](long index)
	{
		return mConcentration[index];
	}
	


	inline const Concentration& ConcentrationSet::operator[](long index) const
	{
		return mConcentration[index];
	}
	
	
	
	inline void ConcentrationSet::SetDeltaFunctionElement(short index, const Element *pElement)
	{
		mDeltaFunctionElement[index] = (Element*) pElement;
		return;
	}
	
	
	
	inline const Element* ConcentrationSet::DeltaFunctionElement(short index) const
	{
		return mDeltaFunctionElement[index];
	}
}


#endif // _concentrationset_h_

