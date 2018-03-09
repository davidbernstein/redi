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

#ifndef _reditypes_h_
#define _reditypes_h_

#include <set>
#include "namespace.h"

namespace NAMESPACE {
	class Reaction;
	class Element;
}

typedef std::list<NAMESPACE::Reaction> ReactionList;

typedef std::set<std::pair<NAMESPACE::Element*, NAMESPACE::Reaction*> > ElementReactionSet;

#endif // _reditypes_h_
