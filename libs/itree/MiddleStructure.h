//------------------------------------------------------------------------------
//  Copyright 2010-2012 by Jyh-Ming Lien and George Mason University
//  See the file "LICENSE" for more information
//------------------------------------------------------------------------------

#pragma once
#ifndef _MIDDLE_STRUCTURE_H_
#define _MIDDLE_STRUCTURE_H_

#include "RectangleTree.h"

template< class rectKD, int K > class MiddleTree;


///////////////////////////////////////////////////////////////////////////////
//
//  MiddleStructure
//
///////////////////////////////////////////////////////////////////////////////

//The middle structure stores the segments containing the node's values
template< class _rectKD, int K >
class MiddleStructure {

public:

    typedef _rectKD rectKD;
    typedef typename rectKD::Int Int;
    typedef typename rectKD::Pt Pt;
    typedef vector<rectKD*> rectList;
    typedef vector<const Pt*> ptList;
    typedef MiddleTree<rectKD,K> Mtree;
    friend class RectangleTree<MiddleStructure,K>; //the only creator of this class

protected: //only RectangleTree can create MiddleStructure
   
    ///////////////////////////////////////////////////////////////////////////
    MiddleStructure(){ 
        msTree_Root=msTree_Left=msTree_Right=NULL;
    }
    
	~MiddleStructure()
	{
		//delete [] msTree_Root;
		if (msTree_Root != msTree_Left) delete msTree_Root;
		delete [] msTree_Left;
	}

    ///////////////////////////////////////////////////////////////////////////
    bool Build( ptList & ptlist, int Level );    
    bool query( rectKD *q, QueryInfo info );
    bool query_lastlevel( rectKD *q, const QueryInfo & info );

    void print(){
        Mtree * leaf=msTree_Left;
        cout<<" * ";
        while( leaf!=NULL ){
            cout<<leaf->getPoint()->getValue()<<" ";
            leaf=leaf->getRight();
        }
        cout<<"\n";
        if( msTree_Root->subTree!=NULL ) msTree_Root->subTree->print();
    }

protected:
    //RectangleTree<rectKD,K> * primary_node;
    //middle sub tree, pointer to root, left most, and right most.
    Mtree * msTree_Root, * msTree_Left, * msTree_Right;
};

///////////////////////////////////////////////////////////////////////////////
//
//  MiddleTree
//
///////////////////////////////////////////////////////////////////////////////
template< class rectKD, int K >
class MiddleTree {

    friend class MiddleStructure<rectKD,K>;

protected:

    typedef typename rectKD::Int Int;
    typedef typename rectKD::Pt Pt;
    typedef MiddleStructure<rectKD,K> MTree;
    typedef typename MTree::rectList rectList;
    //friend class MiddleStructure<rectKD,K>;

    ///////////////////////////////////////////////////////////////////////////
    // Core
    MiddleTree(){
        subTree=NULL;    point=NULL;
        Left=Right=NULL; type=Internal;
		new_right = new_left = false;
    }

	~MiddleTree(){
		delete subTree;
		if (new_right) delete Right;
		if (new_left) delete Left;

		subTree = NULL;
		Left=NULL; 
		Right=NULL;
	}

    bool Build( MiddleTree * leaves, int size, int level );
    MiddleTree * b_search( double pt ); //biniary search on the tree
    bool query(rectKD *q, const QueryInfo & info);
    bool query_sub(rectKD *q, const QueryInfo & info);

    ///////////////////////////////////////////////////////////////////////////
    // Access
    void setAsLeaf(){ type=Leaf; }          //set this node as lead
    void setAsInternal(){ type=Internal; }  //set this node as internal (default)
    void setLeft(MiddleTree *l){ Left=l; }          //set left child
    MiddleTree * getLeft() const { return Left; }   //get left child
    void setRight(MiddleTree *r){ Right=r; }        //set right child
    MiddleTree * getRight() const { return Right; } //get right child
    void setPoint( const Pt * point ){ this->point=point; }
    const Pt * getPoint()const{ return point; }

    ///////////////////////////////////////////////////////////////////////////
    enum NodeType{ Internal=0, Leaf };
    ///////////////////////////////////////////////////////////////////////////
    RectangleTree<MTree,K> * subTree;
    const Pt * point;  //point value for largest value in left subtree
    NodeType type;  //type of this node, either internal or leaf
    MiddleTree * Left;
    MiddleTree * Right;
	bool new_right, new_left;
};

///////////////////////////////////////////////////////////////////////////////
//
//  MiddleStructure_ez
//
///////////////////////////////////////////////////////////////////////////////
template< class rectKD, int K > class MiddleTree_ez;

template< class rectKD, int K >
class MiddleStructure_ez : public MiddleStructure<rectKD,K>
{
    friend class RectangleTree<MiddleStructure_ez,K>; //the only creator of this class

    typedef MiddleStructure<rectKD,K> Base;
    typedef typename Base::Int Int;
    typedef typename Base::Pt Pt;
    typedef typename Base::rectList rectList;
    typedef typename Base::ptList ptList;

protected: //only RectangleTree can create MiddleStructure

	MiddleStructure_ez(){ subTree=NULL; }

    typedef MiddleTree_ez<rectKD,K> Mtree;
    bool Build( ptList & ptlist, int Level );    
    bool query( rectKD *q, QueryInfo info );
	bool query_tree( rectKD *q, const QueryInfo & info );

	void print(){
        Mtree * leaf=(Mtree*)this->msTree_Left;
        cout<<" * ";
        while( leaf!=NULL ){
            cout<<leaf->getPoint()->getValue()<<" ";
            leaf=(Mtree*)(leaf->getRight());
        }
        cout<<"\n";
        if( subTree!=NULL ) subTree->print();
    }

    RectangleTree<MiddleStructure_ez,K> * subTree;
};

///////////////////////////////////////////////////////////////////////////////
//
//  MiddleTree_ez
//
///////////////////////////////////////////////////////////////////////////////
template< class rectKD, int K >
class MiddleTree_ez : public MiddleTree<rectKD,K> {
    
    friend class MiddleStructure_ez<rectKD,K>;

protected:
    typedef MiddleStructure_ez<rectKD,K> MTree;

    bool Build( MiddleTree_ez * leaves, int size, int level );
	bool query(rectKD *q, const QueryInfo & info){return false;}
	bool query_sub(rectKD *q, const QueryInfo & info){return false;}
};

///////////////////////////////////////////////////////////////////////////////
//
//  MiddleStructure Implementation
//
///////////////////////////////////////////////////////////////////////////////
template< class rectKD, int K > bool MiddleStructure<rectKD,K>
::Build( ptList & ptlist, int Level )
{ 
    //check condition
    if( Level>=K || Level<0 ) return false;

    //build middle sub tree.
    int ptSize=ptlist.size();

    Mtree * leaves=new Mtree[ptSize];
    if( leaves==NULL ){ 
        cerr<<"MiddleStructure::Build Error : Not Enough Memory"<<endl; 
        return false; 
    }

    for( int i=0;i<ptSize;i++ ){
        leaves[i].setAsLeaf(); //deafult is internal
        leaves[i].setLeft((i==0)?NULL:&leaves[i-1]);
        leaves[i].setRight((i==ptSize-1)?NULL:&leaves[i+1]);
        leaves[i].setPoint(ptlist[i]);
    }

    msTree_Left=&leaves[0];
    msTree_Right=&leaves[ptSize-1];

    if( ptSize==1 ){ msTree_Root=&leaves[0]; return true; }

    if( (msTree_Root=new Mtree())==NULL ){ 
        cerr<<"MiddleStructure::Build Error : Not Enough Memory"<<endl; 
        return false; 
    }

    if( msTree_Root->Build(leaves,ptSize,Level)==false ){
        cerr<<"MiddleStructure::Build Error"<<endl; 
        return false; 
    }

    return true;
}

template< class rectKD, int K > bool MiddleStructure<rectKD,K>
::query(rectKD *q, QueryInfo info) 
{
    //Mtree * leaf=msTree_Left;
    
    if( info.level==K-1 ){
        return query_lastlevel(q,info);
    }
    else{
        info.level++;
        info.split_found=false;
        return msTree_Root->query(q,info);
    }
    
    return false; //should not happen
}

template< class rectKD, int K > bool MiddleStructure<rectKD,K>
::query_lastlevel(rectKD *q, const QueryInfo & info)
{
    Int & intval=q->getInterval(info.level);
    if ( info.dir == QueryInfo::ALL ){
        Mtree * leaf=msTree_Left;
        while( leaf!=NULL ){
            if (leaf->getPoint()->getType() == Pt::LEFT_ENDPOINT){
                VID id=leaf->getPoint()->getInterval()->getOwner()->getVID();
                q->intersects(id);
            }
            leaf=leaf->getRight();
        }
    }
    else if(info.dir==QueryInfo::LEFT_TO_RIGHT){
        Mtree * leaf=msTree_Root->b_search(intval.getTo().getValue()+1e-10);
        if( (*leaf->getPoint())>intval.getTo() ) leaf=leaf->getLeft();
        while( leaf!=NULL ){
            if (leaf->getPoint()->getType() == Pt::LEFT_ENDPOINT){
                VID id=leaf->getPoint()->getInterval()->getOwner()->getVID();
                q->intersects(id);
            }
            leaf=leaf->getLeft();
        }
    }
    else if(info.dir==QueryInfo::RIGHT_TO_LEFT){
        Mtree * leaf=msTree_Root->b_search(intval.getFrom().getValue()-1e-10);
        if( (*leaf->getPoint())<intval.getFrom() ) leaf=leaf->getRight();
        while( leaf!=NULL ){
            if (leaf->getPoint()->getType() == Pt::RIGHT_ENDPOINT){
                VID id=leaf->getPoint()->getInterval()->getOwner()->getVID();
                q->intersects(id);
            }
            leaf=leaf->getRight();
        }
    }

    return true;
}

///////////////////////////////////////////////////////////////////////////////
//
//  MiddleTree Implementation
//
///////////////////////////////////////////////////////////////////////////////
template< class rectKD, int K > bool MiddleTree<rectKD,K>
::Build( MiddleTree<rectKD,K> * leaves, int size, int level )
{ 
    if( type==Leaf ) return false; //Leaf does not have to "buid"



    int left_size=size/2; int right_size=size-left_size;


    //cout<<"size="<<size<<" MiddleTree<rectKD,K> left_size="<<left_size<<" right_size="<<right_size<<endl;


    if( left_size==1 ){ 
        Left=&leaves[0]; 
        point=Left->getPoint();
    }
    else{ 
        point=leaves[left_size-1].getPoint();
        if( (Left=new MiddleTree())==NULL ) return false;
        Left->Build(leaves,left_size, level);
		new_left = true;
    }

    //cout<<"1 size="<<size<<endl;

    if( right_size==1 ){
        Right=&leaves[left_size];
    }
    else{ 
        if( (Right=new MiddleTree())==NULL ) return false;
        Right->Build(&leaves[left_size],right_size, level);
		new_right = true;
    }

    /*cout<<"A"<<endl;
    cout<<"2 size="<<size<<endl;
    cout<<"K="<<K<<" leve="<<level<<endl;
    cout<<"B"<<endl;
	*/

    //build middle tree
    if( level==K-1 ){
        //cout<<"//there is no more sub tree"<<endl;
        return true; //there is no more sub tree
    }


    rectList rects; rects.reserve(size);
    const Pt * first=leaves[0].getPoint();
    rects.push_back((rectKD*)first->getInterval()->getOwner());
    for( int i=1;i<size;i++ ){
        //cout<<"???@ rects size="<<rects.size()<<" size="<<size<<endl;
        const Pt * point=leaves[i].getPoint();
        if( point->getType()==Pt::RIGHT_ENDPOINT ){
            if( point->getInterval()->getFrom()>=(*first) )
                continue;
        }
        rects.push_back((rectKD*)point->getInterval()->getOwner());
    }

    //cout<<"??? rects size="<<rects.size()<<endl;


    if( (subTree=new RectangleTree<MTree,K>())==NULL ) return false;




    if( subTree->Build(rects,level+1)==false ) return false;

    return true;
}

template< class rectKD, int K > MiddleTree<rectKD,K> * 
MiddleTree<rectKD,K>::b_search( double pt )
{
    if( type==Leaf ) return this;
    if( pt>point->getValue() ) return Right->b_search(pt);
    return Left->b_search(pt);
}

template< class rectKD, int K > bool MiddleTree<rectKD,K>::
query( rectKD *q, const QueryInfo & info )
{
    if( info.dir==QueryInfo::ALL || type==Leaf ) return query_sub(q,info);

    Int & intval=q->getInterval(info.level-1);
    if( info.dir==QueryInfo::LEFT_TO_RIGHT ){
        const Pt & pt=intval.getTo();
        if( (pt.getValue()+1e-10)>point->getValue() ){
            if( Left->query_sub(q,info)==false ) return false; 
            if( Right->query(q,info)==false ) return false;
        }
        else{
            if( Left->query(q,info)==false ) return false;
        }
    }
    else if( info.dir==QueryInfo::RIGHT_TO_LEFT ){
        const Pt & pt=intval.getFrom();
        if( (pt.getValue()-1e-10)<point->getValue() ){
            if( Right->query_sub(q,info)==false ) return false;
            if( Left->query(q,info)==false ) return false;
        }
        else{
            if( Right->query(q,info)==false ) return false;
        }
    }
    else return false; //unknow dir

    return true;
}

template< class rectKD, int K > bool MiddleTree<rectKD,K>::
query_sub(rectKD *q, const QueryInfo & info)
{
    if( type!=Leaf ){
        if( subTree==NULL ) return false;
        return subTree->query(q,info);
    }

    //leaf
    rectKD * rect=(rectKD *)point->getInterval()->getOwner();
    for( int l=info.level-1; l<K; l++ ){
        Int & i1=q->getInterval(l);
        Int & i2=rect->getInterval(l);
        if( i1.getTo()<i2.getFrom() ) return true; //not intersect
        if( i1.getFrom()>i2.getTo() ) return true; //not intersect
    }
    
    q->intersects(rect->getVID());

    return true;
}

///////////////////////////////////////////////////////////////////////////////
//
//  MiddleStructure_ez Implementation
//
///////////////////////////////////////////////////////////////////////////////
template< class rectKD, int K > bool MiddleStructure_ez<rectKD,K>
::Build( ptList & ptlist, int Level )
{ 
    //check condition
    if( Level>=K || Level<0 ) return false;

    //build middle sub tree.
    int ptSize=ptlist.size();

    Mtree * leaves=new Mtree[ptSize];
    if( leaves==NULL ){ 
        cerr<<"MiddleStructure::Build Error : Not Enough Memory"<<endl; 
        return false; 
    }
    for( int i=0;i<ptSize;i++ ){
        leaves[i].setAsLeaf(); //deafult is internal
        leaves[i].setLeft((i==0)?NULL:&leaves[i-1]);
        leaves[i].setRight((i==ptSize-1)?NULL:&leaves[i+1]);
        leaves[i].setPoint(ptlist[i]);
    }

    this->msTree_Left=&leaves[0];
    this->msTree_Right=&leaves[ptSize-1];

    if( ptSize==1 ){ this->msTree_Root=&leaves[0]; return true; }

    if( (this->msTree_Root=new Mtree())==NULL ){
        cerr<<"MiddleStructure::Build Error : Not Enough Memory"<<endl; 
        return false; 
    }

    if( ((Mtree*)this->msTree_Root)->Build(leaves,ptSize,Level)==false ){
        cerr<<"MiddleStructure::Build Error"<<endl; 
        return false; 
    }

	if( Level<K-1 ){
		//build lower dimensional tree
		rectList rects; rects.reserve(ptSize/2);
		typedef typename ptList::iterator IT;
		for(IT iP=ptlist.begin();iP!=ptlist.end();iP++ ){
			if( (*iP)->getType()==Pt::RIGHT_ENDPOINT )  continue;
			rects.push_back((rectKD*)(*iP)->getInterval()->getOwner());
		}
		if( (subTree=new RectangleTree<MiddleStructure_ez,K>())==NULL ) return false;
		if( subTree->Build(rects,Level+1)==false ) return false;
	}
    return true;
}

template< class rectKD, int K > bool MiddleStructure_ez<rectKD,K>
::query(rectKD *q, QueryInfo info) 
{
    query_tree(q,info);

    if( info.level<K-1 ){
        QueryInfo newinfo; newinfo.level=info.level+1;
        if( subTree->query(q,newinfo)==false ) return false;
    }

    return true;
}

template< class rectKD, int K > bool MiddleStructure_ez<rectKD,K>
::query_tree(rectKD *q, const QueryInfo & info)
{ 
	Int & intval=q->getInterval(info.level);
    if ( info.dir == QueryInfo::ALL ){
        Mtree * leaf=(Mtree *)this->msTree_Left;
        while( leaf!=NULL ){
            if (leaf->getPoint()->getType() == Pt::LEFT_ENDPOINT){
                VID id=leaf->getPoint()->getInterval()->getOwner()->getVID();
                q->intersects(id,info.level);
            }
            leaf=(Mtree *)(leaf->getRight());
        }
    }
    else if(info.dir==QueryInfo::LEFT_TO_RIGHT){
        Mtree * leaf=(Mtree *)((Mtree*)this->msTree_Root)->b_search(intval.getTo().getValue()+1e-10);
        if( (*leaf->getPoint())>intval.getTo() ) leaf=(Mtree *)(leaf->getLeft());
        while( leaf!=NULL ){
            if (leaf->getPoint()->getType() == Pt::LEFT_ENDPOINT){
                VID id=leaf->getPoint()->getInterval()->getOwner()->getVID();
                q->intersects(id,info.level);
            }
            leaf=(Mtree *)(leaf->getLeft());
        }
    }
    else if(info.dir==QueryInfo::RIGHT_TO_LEFT){
        Mtree * leaf=(Mtree *)((Mtree*)this->msTree_Root)->b_search(intval.getFrom().getValue()-1e-10);
        if( (*leaf->getPoint())<intval.getFrom() ) leaf=(Mtree *)(leaf->getRight());
        while( leaf!=NULL ){
            if (leaf->getPoint()->getType() == Pt::RIGHT_ENDPOINT){
                VID id=leaf->getPoint()->getInterval()->getOwner()->getVID();
                q->intersects(id,info.level);
            }
            leaf=(Mtree *)(leaf->getRight());
        }
    }

    return true;
}

///////////////////////////////////////////////////////////////////////////////
//
//  MiddleTree_ez Implementation
//
///////////////////////////////////////////////////////////////////////////////

template< class rectKD, int K > bool MiddleTree_ez<rectKD,K>
::Build( MiddleTree_ez<rectKD,K> * leaves, int size, int level )
{ 
    typedef MiddleTree_ez<rectKD,K> me;

    if( this->type==this->Leaf ) return false; //Leaf does not have to "buid"

    int left_size=size/2; int right_size=size-left_size;
    if( left_size==1 ){ 
        this->Left=&leaves[0];
        this->point=((me *)this->Left)->getPoint();
    }
    else{ 
        this->point=leaves[left_size-1].getPoint();
        if( (this->Left=new MiddleTree_ez())==NULL ) return false;
        ((me *)this->Left)->Build(leaves,left_size, level);
    }
    if( right_size==1 ) this->Right=&leaves[left_size];
    else{ 
        if( (this->Right=new MiddleTree_ez())==NULL ) return false;
        ((me *)this->Right)->Build(&leaves[left_size],right_size, level);
    }

    return true;
}

#endif //_MIDDLE_STRUCTURE_H_

