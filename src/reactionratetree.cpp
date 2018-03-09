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

#include "reactionratetree.h"
#include <iostream>

using namespace NAMESPACE;
using namespace std;


Reaction* ReactionRateTree::GetNextReaction(double urn) const
{
	const ReactionRateNode* pNode = Root();
	
	double x = urn * pNode->Value();
	
	
	while (pNode->IsALeaf() == false) {
		pNode = pNode->Search(x);
	}
	
	return pNode->GetReaction();
}


		
void ReactionRateTree::UpdateDependentRates(Reaction *pReaction)
{
	Array<Reaction*>& graph = pReaction->GetDependencyGraph();
	
	for (short i = 0; i < graph.Size(); ++i) 
		Update(graph[i]);
	
	return;
}


		
void ReactionRateTree::Update(const Reaction *pReaction)
{
	ReactionRateNode *pNode = pReaction->GetReactionRateNode();
	
	if (pNode->IsALeaf() == false)
		ThrowException("ReactionRateTree::Update : reaction rate node is not a leaf");
	
	if (pReaction->Rate() == pNode->Value())
		return;
	
	pNode->SetValue(pReaction->Rate());
	pNode = pNode->GetParent();
	while (pNode != NULL) {
		pNode = pNode->SetValueFromChildren();
		//pNode->ComputeValue();
		//pNode = pNode->GetParent();
	}
	
	
	return;
}



void ReactionRateTree::Make(list<Reaction> &reactionList)
{
	// erase current tree
	Erase();
	MakeRoot();
	
	long numReactions = reactionList.size();

	long numLeaves = 1;
	while (2 * numLeaves < numReactions) {
		numLeaves = RefineAllLeaves();
	}
	
	long numLeft = numReactions - numLeaves;
	if (numLeft > numLeaves)
		ThrowException("ReactionRateTree::Make : shouldn't happen");
	
	list<ReactionRateNode>::iterator iNl = mNodeList.begin();
	long i = 0;
	while (i < numReactions - numLeaves) {
		if ((*iNl).IsALeaf()) {
			MakeChildren(&*iNl);
			++i;
		}
		++iNl;
	}
	
	long count = 0;
	for (iNl = mNodeList.begin(); iNl != mNodeList.end(); ++iNl) {
		if ((*iNl).IsALeaf()) {
			++count;
		}
	}
	
	// populate leaves
	list<Reaction>::iterator iRl = reactionList.begin();
	for (iNl = mNodeList.begin(); iNl != mNodeList.end(); ++iNl) {
		if ((*iNl).IsALeaf()) {
			(*iNl).SetReaction(&*iRl);
			(*iRl).SetReactionRateNode(&*iNl);
			++iRl;
		}
	}
	
	SetValues();
	CheckAll();

	return;
}



void ReactionRateTree::SetValues()
{
	SetValues(-1.0);

	list<ReactionRateNode>::iterator iNl;
	for (iNl = mNodeList.begin(); iNl != mNodeList.end(); ++iNl) {
		if ((*iNl).IsALeaf()) {
			Reaction *pReaction = (*iNl).GetReaction();
			(*iNl).SetValue(pReaction->Rate());
		}
	}
	
	// populate rest of the nodes
	ReactionRateNode *pNode = Root();
	while (pNode != NULL) 
		pNode = pNode->ComputeValue();

	
	return;
}



void ReactionRateTree::SetValues(double value)
{
	list<ReactionRateNode>::iterator iNl;
	for (iNl = mNodeList.begin(); iNl != mNodeList.end(); ++iNl)
		(*iNl).SetValue(value);
	
	return;
}



long ReactionRateTree::RefineAllLeaves()
{
	long count = 0;
	list<ReactionRateNode>::iterator iNl = mNodeList.begin();;
	long currentSize = mNodeList.size();
	long i = 0;
	while (i < currentSize) {
		if ((*iNl).IsALeaf()) {
			MakeChildren(&*iNl);
			count += 2;
		}
		++i;
		++iNl;
	}
	
	return count;
}



void ReactionRateTree::MakeChildren(ReactionRateNode *pParent)
{
	ReactionRateNode node;
	mNodeList.push_back(node);
	ReactionRateNode *pChild0 = &*(--mNodeList.end());
	
	mNodeList.push_back(node);
	ReactionRateNode *pChild1 = &*(--mNodeList.end());
	
	pParent->SetChildren(pChild0, pChild1);
	pChild0->SetParent(pParent);
	pChild1->SetParent(pParent);
	
	return;
}



short ReactionRateTree::Depth() const
{	
	short depth = 0;
	list<ReactionRateNode>::const_iterator iNl;
	for (iNl = mNodeList.begin(); iNl != mNodeList.end(); ++iNl) {
		if ((*iNl).IsALeaf()) 
			depth = max(depth, (*iNl).DistanceToRoot());
	}
		
	return depth;
}



void ReactionRateTree::CheckAll() const
{
	list<ReactionRateNode>::const_iterator iNl;
	
	for (iNl = mNodeList.begin(); iNl != mNodeList.end(); ++iNl) {		
		// check reaction pointer
		if ((*iNl).IsALeaf()) {
			Reaction *pReaction = (*iNl).GetReaction();
			if (pReaction == NULL) {
				ThrowException("ReactionRateTree::CheckAll : leaf has NULL reaction pointer");
			}
			else {
				if (&*iNl != pReaction->GetReactionRateNode())
					ThrowException("ReactionRateTree::CheckAll : reaction has wrong ReactionRateNode");
			}
		}
		else {
			if ((*iNl).GetReaction() != NULL) 
				ThrowException("ReactionRateTree::CheckAll : non-leaf has non-NULL reaction pointer");
		}
	
		// check value
		if ((*iNl).IsALeaf() == false) {
			double sum = (*iNl).GetChild(0)->Value() + (*iNl).GetChild(1)->Value();
			double value = (*iNl).Value();
			if (sum != value)
				ThrowException("ReactionRateTree::CheckAll : non-leaf has wrong value");
		}
		else {
			if ((*iNl).Value() < 0.0)
				ThrowException("ReactionRateTree::CheckAll : leaf has negative value");
		}
		
		// check element
		if ((*iNl).IsALeaf()) {
			Element *pE = (*iNl).GetReaction()->GetElement();
			pE->GetOccupancyData();
		}
	}
	
	
	return;
}



void ReactionRateTree::Print() const
{
	cout << "ReactionRateTree Print start" << endl;
	
	list<ReactionRateNode>::const_iterator iNl;
	for (iNl = mNodeList.begin(); iNl != mNodeList.end(); ++iNl) {
		(*iNl).Print();
	}
	
	cout << "ReactionRateTree Print end" << endl;


	return;
}
