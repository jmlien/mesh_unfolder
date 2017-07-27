//------------------------------------------------------------------------------
//  Copyright 2010-2012 by Jyh-Ming Lien and George Mason University
//  See the file "LICENSE" for more information
//------------------------------------------------------------------------------

#ifndef _ITREE_RECKD_H_
#define _ITREE_RECKD_H_

#ifdef WIN32
#pragma warning( disable : 4786 )
#endif

///////////////////////////////////////////////////////
// STL Headers
#include <math.h>
#include <string>
#include <string.h>
#include <list>
#include <vector>
#include <algorithm>
#include <iostream>
using namespace std;

///////////////////////////////////////////////////////
typedef unsigned int VID;
template<class Pt>  class Interval;

///////////////////////////////////////////////////////////////////////////////
//
//  Rect
//
///////////////////////////////////////////////////////////////////////////////

class Rect 
{

public:

    typedef vector< VID > Intersect;

    Rect(VID _vid)
    {
        vid = _vid;
        isize=0;
        magic_byte=0;
    }

    ///////////////////////////////////////////////////////
    // Core
    void reportIntersections();

    ///////////////////////////////////////////////////////
    // ACCESS
    void setVID(VID vid) { this->vid = vid; }
    const VID getVID() const { return vid; }
    Intersect & getIntersections() { return intersection; }
    const Intersect & getIntersections() const { return intersection; }

    char magic_byte;
    
protected:

    //static unsigned int rect_size; //record how many of the rect is created.
    
    Intersect intersection;
    VID vid;
    
    int isize;//temp use, should delete later
};

///////////////////////////////////////////////////////////////////////////////
//
//  RectKD
//
///////////////////////////////////////////////////////////////////////////////

template < class Point, int K >
class RectKD : public Rect  //a K dimensional rectangle
{ 
public:

	typedef Point Pt;
    typedef Interval<Pt> Int;

    ///////////////////////////////////////////////////////
    // Con / Des
	RectKD(VID id) : Rect(id){}
    void intersects( VID id ){
        if( id==vid ) return;
        isize++;
        intersection.push_back(id);
    }

    ///////////////////////////////////////////////////////
    // ACCESS
    void setInterval( const Int & intv, int k ){ 
        if( k<0 ) k=0;
        if( k>=K ) k=K-1;
        m_Intervals[k]=intv;
        m_Intervals[k].setOwner(this);
    }

    Int & getInterval( int k ){ 
        if( k<0 ) k=0;
        if( k>=K ) k=K-1;
        return m_Intervals[k];
    }

private:

    Int m_Intervals[K];
};

///////////////////////////////////////////////////////////////////////////////
//
//  RectKD_ez
//
///////////////////////////////////////////////////////////////////////////////
template < class Pt, int K > //K is dimension
class RectKD_ez : public RectKD<Pt,K>   //a K dimensional rectangle
{

public:

    typedef short BYTE;

    ///////////////////////////////////////////////////////
    // Con / Des
	RectKD_ez(){ m_Result=NULL; }
	~RectKD_ez(){ delete[] m_Result; }

    bool beginQuery(){
        static unsigned int cellsize=(int)ceil( (double)(K*this->rect_size)/(8*sizeof(BYTE)) );
        if( (m_Result=new BYTE[cellsize])==NULL ) return false;
		memset(m_Result,0,cellsize*sizeof(BYTE));
        return true;
    }

    void intersects( VID id, int level ){
        static unsigned int byte_size=sizeof(BYTE)*8;
		if( m_Result==NULL ){ 
			if( beginQuery()==false ){
				cerr<<"RectKD_ez Error:: Create m_Result Error"<<endl;
				exit(1);
			}
		}
		int index=(id*K)+level;
        int byte_id=index/byte_size; 
        int bit_id=index-byte_size*byte_id;
        m_Result[byte_id]=m_Result[byte_id] | (1<<bit_id);
    }

	void reportIntersections(){
		endQuery();
		RectKD<Pt,K>::reportIntersections();
	}

    bool endQuery(){
        if( m_Result==NULL ) return false;
        static unsigned int byte_size=sizeof(BYTE)*8;
        for( unsigned int iR=0;iR<this->rect_size;iR++ ){
            if( iR==this->vid ) continue;
			bool intersect=true;
			int index=iR*K;

			for( int iK=0;iK<K;iK++ ){
				int id=index+iK;
			    int byte_id=id/byte_size; 
		        int bit_id=id-byte_size*byte_id;
				if( (m_Result[byte_id] & (1<<bit_id))==0 ){
					intersect=false;
					break;
				}
			}
            if( intersect ) this->intersection.push_back(iR);
        }
        delete [] m_Result;
        m_Result=NULL;
        return true;
    }

private:
    BYTE * m_Result;
};


template < class Pt, int K >
ostream & operator << (ostream & out, const RectKD<Pt,K>& rectKD);


#endif //_ITREE_RECKD_H_
