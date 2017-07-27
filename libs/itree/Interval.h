//------------------------------------------------------------------------------
//  Copyright 2010-2012 by Jyh-Ming Lien and George Mason University
//  See the file "LICENSE" for more information
//------------------------------------------------------------------------------

#pragma once
#ifndef _ITREE_INTERVAL_H_
#define _ITREE_INTERVAL_H_

///////////////////////////////////////////////////////////////////////////////
#include <set>
#include <iostream>
using namespace std;

///////////////////////////////////////////////////////////////////////////////
class Rect;

template<class Pt>
class Interval {
public:

    ///////////////////////////////////////////////////////////////////////////////
    //constructors:
    Interval(){ owner=NULL; } //do nothing
    Interval( const Interval<Pt> & other ){ *this=other; }
    Interval( Rect * owner,const Pt& _from, const Pt& _to)
    {   //create endpts
        this->owner=owner;
        setFrom(_from);
        setTo(_to);
    }
    
    ///////////////////////////////////////////////////////////////////////////////
    //Access
    const Pt & getFrom() const { return from; }
    void setFrom(const Pt& _from) { 
        from = _from; 
        from.setOwner(this);
        from.setType(Pt::LEFT_ENDPOINT);
    }
    const Pt & getTo() const { return to; }
    void setTo(const Pt& _to) { 
        to = _to; 
        to.setOwner(this);
        to.setType(Pt::RIGHT_ENDPOINT);
    }

	Rect * getOwner() { return owner; }
    void setOwner(Rect * owner) { this->owner=owner; }

    ///////////////////////////////////////////////////////////////////////////////
    //Core
    enum positions {RIGHT, LEFT, IN, UNKNOWN};

    //find position of given endpt related to this interval
    enum positions position(const Pt& value){
        if (from < value && to < value) return RIGHT;
        else if (from > value && to > value) return LEFT;
        else if (from <= value && to >= value) return IN;
        else return UNKNOWN;
    }

    const Interval<Pt> & operator=( const Interval<Pt> & other ) {
        setFrom(other.from); setTo(other.to);
        owner=other.owner;
        return *this;
    }
    
private:
    Pt from; //origin of the interval
    Pt to;   //span of the interval
    Rect * owner; //polyhedron which the interval belongs to
};

///////////////////////////////////////////////////////////////////////////////
template< class Pt >
ostream & operator << (ostream & out, const Interval<Pt> & interval);

#endif //_ITREE_INTERVAL_H_
