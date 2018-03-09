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

#ifndef _reaction_h_
#define _reaction_h_

#include "element.h"
#include "reditypes.h"

namespace NAMESPACE {
	class ReactionRateNode;
}

namespace NAMESPACE {
	class Reaction {	
	public:
		// Constructor
		Reaction(void);

		// Copy constructor
		Reaction(const Reaction &r);

		// Destructor
		~Reaction(void);

		// operators
		Reaction& operator=(const Reaction &r);
		
		// reaction index
		void SetReactionIndex(short value);
		short ReactionIndex(void) const;
		void SetNeighborIndex(short value);
		short NeighborIndex(void) const;
		
		// rate
		void SetRate(double value);
		double Rate(void) const;
		
		// element
		void SetElement(const Element* pE);
		Element* GetElement(void) const;
		
		// dependency graph
		void InitializeDependencyGraph(std::set<Reaction*> &reactionSet);
		bool IsADependentReaction(const Reaction *pReaction) const;
		short NumDependentReactions(void) const;
		Array<Reaction*>& GetDependencyGraph(void);
		
		// reaction-rate node
		void SetReactionRateNode(const ReactionRateNode* pReactionRateNode);
		ReactionRateNode* GetReactionRateNode(void) const;
		
		// set iterator
		void SetIterator(void *pReactionMultiSetIterator);
		void* GetIterator(void) const;

		// IO
		void Print(void) const;
		
		// operators needed for STL
		bool operator==(const Reaction &r) const;
		bool operator<(const Reaction &r) const;
		bool operator!=(const Reaction &r) const;
		bool operator>(const Reaction &r) const;
		
		// private member functions
	private:
		

	private:
		// reaction index (which reaction in chemicalSysatem this is)
		short mReactionIndex[2];
		
		// rate
		double mRate;
		
		// element
		Element *mpElement;
		
		// dependent reactions
		Array<Reaction*> mDependentReaction;
		
		// ReactionRateNode
		ReactionRateNode* mpReactionRateNode;
		
		// position on list
		void *mpListIterator;
	};



	inline Reaction::Reaction()
	{
		mpListIterator = NULL;
		mRate = -1.0;
		mReactionIndex[0] = -1;
		mReactionIndex[1] = -1;
		mpReactionRateNode = NULL;
		
		return;
	}
	
	
	
	inline Reaction::~Reaction()
	{
    	if (mpListIterator != NULL) {
        	delete ((std::list<Reaction>::iterator *) mpListIterator);
            mpListIterator = NULL;
    	}
		
        return;
	}
	
	
	
	inline Reaction::Reaction(const Reaction &r)
	{
		// copy everything but iterator		
		mReactionIndex[0] = r.mReactionIndex[0];
		mReactionIndex[1] = r.mReactionIndex[1];
		mRate = r.mRate;
		mpElement = r.mpElement;
		mDependentReaction = r.mDependentReaction;
		
		mpListIterator = NULL;
		
		return;
	}
	
	
	
	inline void Reaction::SetReactionIndex(short value)
	{
		mReactionIndex[0] = value;
		return;
	}
	
	
	
	inline short Reaction::ReactionIndex() const
	{
		return mReactionIndex[0];
	}
	
	
	
	inline void Reaction::SetNeighborIndex(short value)
	{
		mReactionIndex[1] = value;
		return;
	}
	
	
	
	inline short Reaction::NeighborIndex() const
	{
		return mReactionIndex[1];
	}
	
	
	
	inline void Reaction::SetRate(double value)
	{
		mRate = value;
		return;
	}
	
	
	
	inline double Reaction::Rate() const
	{
		return mRate;
	}	
	
	
	
	inline void Reaction::SetElement(const Element* pE)
	{
		mpElement = (Element*) pE;
		return;
	}
	
	
	
	inline Element* Reaction::GetElement() const
	{
		return mpElement;
	}
	
	
	
	inline Array<Reaction*>& Reaction::GetDependencyGraph()
	{
		return mDependentReaction;
	}
	
	
	
	inline void Reaction::InitializeDependencyGraph(std::set<Reaction*> &reactionSet)
	{
		mDependentReaction = reactionSet;
		return;
	}
	
	
	
	inline bool Reaction::IsADependentReaction(const Reaction *pReaction) const
	{
		short i = 0;
		while (i < mDependentReaction.Size()) {
			if (pReaction == mDependentReaction[i])
				return true;
			
			++i;
		}
		
		return false;
	}
	
	
	
	inline short Reaction::NumDependentReactions() const
	{
		return mDependentReaction.Size();
	}



	inline void* Reaction::GetIterator() const
	{
        return mpListIterator;
	}
	
	
	
	inline void Reaction::SetIterator(void *pIterator)
	{
		if (mpListIterator != NULL)
			ThrowException("Reaction::SetIterator : already active");
			
		mpListIterator = new ReactionList::iterator;
		*(ReactionList::iterator *) mpListIterator = *(ReactionList::iterator *) pIterator;

		return;
	}
	
	
	
	inline void Reaction::SetReactionRateNode(const ReactionRateNode* pReactionRateNode)
	{
		mpReactionRateNode = (ReactionRateNode*) pReactionRateNode;
		return;
	}
	
	
	
	inline ReactionRateNode* Reaction::GetReactionRateNode() const
	{
		return mpReactionRateNode;
	}
}


#endif // _reaction_h_

