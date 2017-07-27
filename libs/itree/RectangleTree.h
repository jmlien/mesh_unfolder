//------------------------------------------------------------------------------
//  Copyright 2010-2012 by Jyh-Ming Lien and George Mason University
//  See the file "LICENSE" for more information
//------------------------------------------------------------------------------

#pragma
#ifndef _RECTANGLETREE_H_
#define _RECTANGLETREE_H_

///////////////////////////////////////////////////////////
// Include stl headers
#include <vector>
#include <iostream>
#include <algorithm>
#include <list>
using namespace std;

///////////////////////////////////////////////////////////
// Include itree headers
#include "Interval.h"
#include "EndPoint.h"
#include "RectKD.h"

///////////////////////////////////////////////////////////////////////////////
//
//  QueryInfo
//
///////////////////////////////////////////////////////////////////////////////

//data for performing query
class QueryInfo{
public:
    enum Branch {LEFT_BRANCH, RIGHT_BRANCH};
    enum query_direction{LEFT_TO_RIGHT, RIGHT_TO_LEFT, ALL};

    QueryInfo(){ 
        level=0; split_found=false; 
        branch=LEFT_BRANCH; dir=LEFT_TO_RIGHT;
    }

    int level;
    bool split_found;
    Branch branch;
    query_direction dir;
};

///////////////////////////////////////////////////////////////////////////////
//
//  RectangleTree
//
///////////////////////////////////////////////////////////////////////////////

//The rectangle tree is the data structure storing intervals to be used in 
//queries for intersections
template< class MTree, int K >
class RectangleTree {

    typedef typename MTree::Int Int;
    typedef typename MTree::Pt Pt;
    typedef typename MTree::rectKD rectKD;
    typedef typename MTree::rectList rectList;
    typedef typename MTree::ptList ptList;

public:
    
    ///////////////////////////////////////////////////////////////////////////
    RectangleTree(){ left=right=NULL; mTree=NULL; value=0.0; }
	~RectangleTree()
	{ 
		delete left; 
		delete right; 
		delete mTree; 
	}
    
    ///////////////////////////////////////////////////////////////////////////
    bool Build( rectList& rects, int Level=0 );
    bool query( rectKD *q ){ return query(q,QueryInfo()); }
    bool query( rectKD *q, QueryInfo info );

    void print(){
        if( mTree!=NULL ) mTree->print();
        if( left!=NULL ) left->print();
        if( right!=NULL ) right->print();
    }

///////////////////////////////////////////////////////////////////////////
//Protected and Private
protected:

    bool BuildTree( ptList & ptlist, int Level );
    bool BuildLeft( ptList & ptlist, int Level );
    bool BuildMiddle( ptList & ptlist, int Level );
    bool BuildRight( ptList & ptlist, int Level );

    bool query_no_split( rectKD *q, QueryInfo & info );
    bool query_have_split( rectKD *q, QueryInfo & info );
    bool queryAllSub( rectKD* q, QueryInfo info );

private:
    double value;
    RectangleTree<MTree,K> *left;
    RectangleTree<MTree,K> *right;
    MTree * mTree;
};

///////////////////////////////////////////////////////////////////////////////
//
//  RectangleTree Implementation
//
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
//
// Build R Tree
//
///////////////////////////////////////////////////////////////////////////////

template< class MTree, int K > bool RectangleTree<MTree,K>
::Build( rectList& rects, int Level )
{
    if( Level>=K || Level<0 ) return false; //wrong level

    int size=rects.size();

    ptList endPt; endPt.reserve(size*2);

    //get all end points and sort them
    for(typename rectList::iterator i=rects.begin(); i!=rects.end(); ++i )
    {
        Int & interval= (*i)->getInterval(Level);
        endPt.push_back(&interval.getFrom());
        endPt.push_back(&interval.getTo());
    }

    sort(endPt.begin(), endPt.end(), EndPoint_less());

    return BuildTree(endPt, Level);
}

template< class MTree, int K > bool RectangleTree<MTree,K>
::BuildTree( ptList & ptlist, int Level )
{ 
    //check condition
    if( Level>=K || Level<0 ) return false;
    unsigned int pt_size = ptlist.size();
    if( pt_size<2 ) return false;       //this should not happen (ptlist should have at least 2)
    if( pt_size%2==1) return false;     //this should not happen (ptlist should have even size)
    
    //compute Sl, Sr, and Sm
    ptList Ptleft; ptList Ptright; ptList Ptmiddle;

    //find mid value
    value=ptlist[pt_size/2]->getValue();

    //distinguish left, right , middle
    typedef typename ptList::iterator IT;

    int left_count=0;
    for(IT iP=ptlist.begin();iP!=ptlist.end();++iP){
        if( (*iP)->getType()==Pt::LEFT_ENDPOINT ) left_count++;
    }

    for(IT iP=ptlist.begin();iP!=ptlist.end();++iP)
    {
        Int * interval= (*iP)->getInterval();
        rectKD * rect=(rectKD *)interval->getOwner();
        typename Int::positions pos;

        if( (*iP)->getType()==Pt::LEFT_ENDPOINT ){
            pos=interval->position(value);
            rect->magic_byte=(char)pos;
        }
        else{
            pos=(typename Int::positions)rect->magic_byte;
        }

        switch(pos) {
            //value is at right of i, so put i in Rleft
            case Int::LEFT:   Ptright.push_back(*iP);  break;
            case Int::RIGHT:  Ptleft.push_back(*iP); break;
            case Int::IN:     Ptmiddle.push_back(*iP); break;
            case Int::UNKNOWN: assert(false); break;
        }
    }

    ptlist.clear(); //dose this reduce the memory required?

    //build children
    if( BuildLeft( Ptleft, Level )==false ) return false;
    if( BuildMiddle( Ptmiddle, Level )==false ) return false;
    if( BuildRight( Ptright, Level )==false ) return false;

    return true;
}

template< class MTree, int K > bool RectangleTree<MTree,K>::
BuildLeft( ptList & ptlist, int Level )
{
    if (!ptlist.empty())
    {
        if( (left=new RectangleTree<MTree,K>())==NULL ) return false;
        if( left->BuildTree(ptlist,Level)==false ) return false;
    }

    return true;
}

template< class MTree, int K > bool RectangleTree<MTree,K>::
BuildRight( ptList & ptlist, int Level )
{
    if (!ptlist.empty())
    {
        if( (right=new RectangleTree<MTree,K>())==NULL ) return false;
        if( right->BuildTree(ptlist,Level)==false ) return false;
    }

    return true;
}

template< class MTree, int K > bool RectangleTree<MTree,K>::
BuildMiddle( ptList & ptlist, int Level )
{
    if (!ptlist.empty())
    {
        if( (mTree=new MTree())==NULL ) return false;
        if( mTree->Build(ptlist,Level)==false ) return false;
    }

    return true;
}

///////////////////////////////////////////////////////////////////////////////
//
// Query
//
///////////////////////////////////////////////////////////////////////////////
template< class MTree, int K > bool RectangleTree<MTree,K>::
query( rectKD *q, QueryInfo info )
{
    if (!info.split_found)  return query_no_split(q, info);
    return query_have_split(q, info);
}

template< class MTree, int K > bool RectangleTree<MTree,K>::
query_no_split( rectKD *q, QueryInfo & info )
{
    //find value's position relative to q
    Int & intval=q->getInterval(info.level);
    typename Int::positions pos=intval.position(value);

    switch(pos) {
    case Int::RIGHT:

        if( mTree!=NULL ){
			info.dir=QueryInfo::LEFT_TO_RIGHT;
            if( mTree->query(q, info)==false ) return false;
		}
        if( left!=NULL ) if( left->query(q, info)==false ) return false;
        break;
    case Int::LEFT: 

        if( mTree!=NULL ){
			info.dir=QueryInfo::RIGHT_TO_LEFT;
            if( mTree->query(q, info)==false ) return false;
		}
        if( right!=NULL ) right->query(q, info);
        break;
    case Int::IN: default:

        if( mTree!=NULL ){
			info.dir=QueryInfo::ALL;
            if( mTree->query(q, info)==false ) return false;
		}
        info.split_found=true;
        if( left!=NULL ){
            info.branch=QueryInfo::LEFT_BRANCH;
            if( left->query(q, info)==false ) return false;
        }
        if( right!=NULL ){
            info.branch=QueryInfo::RIGHT_BRANCH;
            if( right->query(q, info)==false ) return false;
        }
        break;
    }

    return true;
}

template< class MTree, int K > bool RectangleTree<MTree,K>::
query_have_split( rectKD *q, QueryInfo & info )
{
    //find value's position relative to q
    Int & intval=q->getInterval(info.level);
    typename Int::positions pos=intval.position(value);

    if (info.branch == QueryInfo::RIGHT_BRANCH) {
        switch(pos) {
        case Int::RIGHT:
            if( mTree!=NULL ){
				info.dir=QueryInfo::LEFT_TO_RIGHT;
                if( mTree->query(q, info)==false ) return false;
			}
            if( left!=NULL ) if( left->query(q, info)==false ) return false;
            break;
        case Int::IN: default:
            if( mTree!=NULL ) {
				info.dir=QueryInfo::ALL;
                if( mTree->query(q, info)==false ) return false;
			}
            if( left!=NULL ) if( left->queryAllSub(q, info)==false ) return false;
            if( right!=NULL ) if( right->query(q, info)==false ) return false;
            break;
        }
    }
    else { //branch == QueryInfo::EFT_BRANCH
        switch(pos) {
        case Int::LEFT:
            if( mTree!=NULL ){
				info.dir=QueryInfo::RIGHT_TO_LEFT;
                if( mTree->query(q, info)==false ) return false;
			}
            if( right!=NULL ) if( right->query(q, info)==false ) return false;
            break;
        case Int::IN: default:
            if( mTree!=NULL ){
				info.dir=QueryInfo::ALL;
                if( mTree->query(q, info)==false ) return false;
			}
            if( right!=NULL ) if( right->queryAllSub(q,info)==false ) return false;
            if( left!=NULL ) if( left->query(q,info)==false ) return false;
            break;
        }
    }

    return true;
}

template< class MTree, int K > bool RectangleTree<MTree,K>
::queryAllSub( rectKD* q, QueryInfo info )
{
    if( left!=NULL ) if( left->queryAllSub(q,info)==false ) return false;
    if( mTree!=NULL ){
		info.dir=QueryInfo::ALL;
        if( mTree->query(q, info)==false ) return false;
	}
    if( right!=NULL ) if( right->queryAllSub(q,info)==false ) return false;

    return true;
}

#endif //_RECTANGLETREE_H_
