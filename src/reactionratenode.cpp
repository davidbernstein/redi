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

#include "reactionratenode.h"
#include <iostream>

using namespace NAMESPACE;
using namespace std;

	
void ReactionRateNode::Print() const
{
	cout << "this = " << this << endl;
	cout << "parent = " << mpParent << endl;
	
	if (mpChild != NULL) {
		cout << "child 0 = " << &mpChild[0] << endl;
		cout << "child 1 = " << &mpChild[1] << endl;
	}
	
	cout << "value = " << mValue << endl;
	cout << "reaction pointer " << mpReaction << endl;
	
	cout << endl;
	
	return;
}		



short ReactionRateNode::DistanceToRoot() const
{
	short d = 0;
	const ReactionRateNode *pNode = mpParent;
	while (pNode != NULL) {
		++d;
		pNode = pNode->mpParent;
	}
	
	return d;
}

