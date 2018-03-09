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

#ifndef _reactionratenode_h_
#define _reactionratenode_h_

#include "reaction.h"

namespace NAMESPACE {
	class ReactionRateNode {	
	public:
		// Constructor
		ReactionRateNode(void);

		// Copy constructor
		ReactionRateNode(const ReactionRateNode &n);

		// Destructor
		~ReactionRateNode(void);
		void Erase(void);
		
		// parents and children
		ReactionRateNode* GetParent(void) const;
		void SetParent(ReactionRateNode *pParent);
		void SetChildren(const ReactionRateNode *pChild0, const ReactionRateNode *pChild1);
		bool IsALeaf(void) const;
		ReactionRateNode* GetChild(short index) const;
		ReactionRateNode* LeftChild(void) const;
		ReactionRateNode* RightChild(void) const;
		
		// value
		void SetValue(double value);
		double Value(void) const;
		ReactionRateNode* ComputeValue(void);
		ReactionRateNode* SetValueFromChildren(void);
		
		// reaction
		void SetReaction(const Reaction *pReaction);
		Reaction* GetReaction(void) const;
		
		// searching
		ReactionRateNode* Search(double &x) const;
		short DistanceToRoot(void) const;
		
		// IO
		void Print(void) const;

	private:
		// pointer to its parent
		ReactionRateNode *mpParent;

		// children
		ReactionRateNode **mpChild;
		
		// value
		double mValue;
		
		// pointer to corresponding reaction
		Reaction *mpReaction;
	};



	inline ReactionRateNode::ReactionRateNode()
	{
		mpParent = NULL;
		mpChild = NULL;
		mValue = -1.0;
		mpReaction = NULL;

		return;
	}
	
	
	
	inline ReactionRateNode::ReactionRateNode(const ReactionRateNode &n)
	{
		mpParent = n.mpParent;
		mValue = n.mValue;
		mpReaction = n.mpReaction;
		
		if (n.mpChild == NULL) {
			mpChild = NULL;
		}
		else {
			if (mpChild == NULL)
				mpChild = new ReactionRateNode*[2];
			
			mpChild[0] = n.mpChild[0];
			mpChild[1] = n.mpChild[1];
		}


		return;
	}
	


	inline ReactionRateNode::~ReactionRateNode()
	{
		Erase();
		return;
	}
	
	
	inline void ReactionRateNode::Erase()
	{
		if (mpChild != NULL) {
			delete [] mpChild;
			mpChild = NULL;
		}
		
		mpParent = NULL;
		mValue = -1.0;
		mpReaction = NULL;
		
		return;
	}
	


	inline void ReactionRateNode::SetParent(ReactionRateNode *pParent)
	{
		mpParent = pParent;
		return;
	}



	inline ReactionRateNode* ReactionRateNode::GetParent() const
	{
		return mpParent;
	}
	
	
	
	inline bool ReactionRateNode::IsALeaf() const
	{
		return mpChild == NULL;
	}
	
	
	
	inline void ReactionRateNode::SetValue(double value)
	{
		mValue = value;
		return;
	}
	
	
	
	inline double ReactionRateNode::Value(void) const
	{
		return mValue;
	}
	
	
	
	inline ReactionRateNode* ReactionRateNode::ComputeValue()
	{
		if (mpChild != NULL) {
			if (mpChild[0]->Value() < 0.0)
				return mpChild[0];
			
			if (mpChild[1]->Value() < 0.0)
				return mpChild[1];
			
			mValue = mpChild[0]->Value() + mpChild[1]->Value();
			
			return mpParent;
		}
		else {
			// make sure my value has been set
			if (mValue < 0.0)
				ThrowException("ReactionRateNode::ComputeValue : leaf value not set");
			
			return mpParent;
		}
	}
	
	
	
	inline ReactionRateNode* ReactionRateNode::SetValueFromChildren()
	{
		if (mpChild == NULL)
			ThrowException("ReactionRateNode::SetValueFromChildren : Node is a leaf");
		
		mValue = mpChild[0]->Value() + mpChild[1]->Value();
		
		return mpParent;
	}
	 
		
		
	inline void ReactionRateNode::SetReaction(const Reaction *pReaction)
	{
		mpReaction = (Reaction*) pReaction;
		return;
	}
	
	
	
	inline Reaction* ReactionRateNode::GetReaction(void) const
	{
		return mpReaction;
	}
	
	
	
	inline void ReactionRateNode::SetChildren(const ReactionRateNode *pChild0, const ReactionRateNode *pChild1)
	{
		if (mpChild == NULL) 
			mpChild = new ReactionRateNode*[2];
		
		mpChild[0] = (ReactionRateNode*) pChild0;
		mpChild[1] = (ReactionRateNode*) pChild1;
			
		mpChild[0]->SetParent(this);
		mpChild[1]->SetParent(this);
		
		return;
	}
	
	
	
	inline ReactionRateNode* ReactionRateNode::GetChild(short index) const
	{
		if (mpChild == NULL)
			return NULL;
		
		return mpChild[index];
	}
		
		
		
	inline ReactionRateNode* ReactionRateNode::Search(double &x) const
	{
		if (x < mpChild[0]->Value()) {
			return mpChild[0];
		}
		else {
			x -= mpChild[0]->Value();
			return mpChild[1];
		}
	}

	
	
	inline ReactionRateNode* ReactionRateNode::LeftChild() const
	{
		return mpChild[0];
	}
	
	
	
	inline ReactionRateNode* ReactionRateNode::RightChild() const
	{
		return mpChild[1];
	}
}


#endif // _reactionratenode_h_

