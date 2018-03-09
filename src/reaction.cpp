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

#include "reaction.h"

#include <iostream>

using namespace NAMESPACE;
using namespace std;

void Reaction::Print() const
{
	cout << "this = " << this << endl;
	cout << "value = " << mRate << endl;
	cout << "reaction index = " << mReactionIndex[0] << endl;
	cout << "neighbor index = " << mReactionIndex[1] << endl;
	cout << "element pointer = " << mpElement << endl;
	
	if (mReactionIndex[1] >= 0)
		cout << "neighbor element pointer = " << mpElement->PN(mReactionIndex[1]) << endl;
	
	cout << "reaction rate node pointer = " << endl;
	cout << endl;

	return;
}



bool Reaction::operator==(const Reaction &r) const
{
	ThrowException("Reaction::operator== called");
	return false;
}



bool Reaction::operator!=(const Reaction &r) const
{
	ThrowException("Reaction::operator!= called");
	return false;
}



bool Reaction::operator<(const Reaction &r) const
{
	ThrowException("Reaction::operator< called");
	return false;
}



bool Reaction::operator>(const Reaction &r) const
{
	ThrowException("Reaction::operator> called");
	return false;
}
