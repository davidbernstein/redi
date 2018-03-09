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

#ifndef _reactionratetree_h_
#define _reactionratetree_h_

#include "reactionratenode.h"
#include <list>

namespace NAMESPACE {
	class ReactionRateTree {	
	public:
		// Constructor
		ReactionRateTree(void) { };

		// Destructor
		~ReactionRateTree(void);
		void Erase(void);
			
		// initialization
		void MakeRoot(void);
		void Make(std::list<Reaction> &reactionList);
		void MakeChildren(ReactionRateNode *pParent);
		void SetValues(void);
		void SetValues(double value);
		
		// update
		void UpdateDependentRates(Reaction *pReaction);
		void Update(const Reaction *pReaction);
		
		// find
		Reaction* GetNextReaction(double urn) const;
		double GetNextReactionTime(double urn) const;
	
		// consistency checking
		void CheckAll(void) const;
		void Print(void) const;
		short Depth(void) const;
		ReactionRateNode* Root(void) const;
		
	private:
		//void MakeLeafList(std::list<ReactionRateNode*> &leafList) const;
		long RefineAllLeaves(void);

	private:
		std::list<ReactionRateNode> mNodeList;
	};
	
	
	
	inline ReactionRateTree::~ReactionRateTree()
	{
		mNodeList.erase(mNodeList.begin(), mNodeList.end());
		return;
	}
	
	
	
	inline void ReactionRateTree::Erase()
	{
		mNodeList.erase(mNodeList.begin(), mNodeList.end());
		return;
	}
	
	
	
	inline void ReactionRateTree::MakeRoot()
	{
		ReactionRateNode node;
		mNodeList.push_back(node);
	}
	
	
	
	inline double ReactionRateTree::GetNextReactionTime(double urn) const
	{
		return -log(urn) / Root()->Value();
	}
	
	
	
	inline ReactionRateNode* ReactionRateTree::Root() const 
	{
		return (ReactionRateNode*) &*mNodeList.begin();
	}
}


#endif // _reactionratetree_h_

