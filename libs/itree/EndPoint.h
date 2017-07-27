//------------------------------------------------------------------------------
//  Copyright 2010-2012 by Jyh-Ming Lien and George Mason University
//  See the file "LICENSE" for more information
//------------------------------------------------------------------------------


//The endpoints are used in the middle structures to relate the endpoints of the 
//intervals with the intervals and the polyhedra they belong to.

#pragma once
#ifndef _ITREE_ENDPOINT_H_
#define _ITREE_ENDPOINT_H_

///////////////////////////////////////////////////////////////////////////////
#include <iostream>
using namespace std;

///////////////////////////////////////////////////////////////////////////////
template< class Pt > class Interval;

///////////////////////////////////////////////////////////////////////////////
// End Point of Interval. double need have operator < and > for compare and 
// << for output and = for assignment
class EndPoint {
public:
    enum endpointType {LEFT_ENDPOINT, RIGHT_ENDPOINT};
    typedef Interval< EndPoint > Int;

    ///////////////////////////////////////////////////////////////////////////////
    // Con / Des
    EndPoint()
    { 
        owner=NULL; endpoint_type=RIGHT_ENDPOINT; 
    }

    EndPoint( const double & value )
    { 
        setValue(value);
        owner=NULL;
        endpoint_type=RIGHT_ENDPOINT;
    }

    EndPoint(Int * owner, const double & value, endpointType endpoint_type=RIGHT_ENDPOINT)
    {
        setOwner(owner); setType(endpoint_type); setValue(value);
    }

    EndPoint(const EndPoint& pt){ *this=pt; }

    ///////////////////////////////////////////////////////////////////////////////
    // Access
    void setValue( const double & v ) { value=v; }
    double getValue() const { return value; }

    void setOwner(Int * owner){ 
        this->owner=owner; 
    }
    Int * getInterval() const { return owner; }

    void setType(endpointType type) { endpoint_type=type; }
    endpointType getType() const { return endpoint_type; }

    ///////////////////////////////////////////////////////////////////////////////
    // Operator
    bool operator < (const EndPoint& other) const{return other.value>value;}
    bool operator <=(const EndPoint& other) const{ return !(*this>other); }
    bool operator > (const EndPoint& other) const{return value>other.value;}
    bool operator >=(const EndPoint& other) const{ return !(*this<other); }
    bool operator ==(const EndPoint& other) const{ return value==other.value; }
    EndPoint operator+(const EndPoint& other) const{ return EndPoint(value+other.value); }
    EndPoint operator/(const EndPoint& other ) const{ return EndPoint(value/other.value); }
    const EndPoint & operator = (const EndPoint& other) { 
        value=other.value; owner=other.owner; endpoint_type=other.endpoint_type;
        return *this;
    }

///////////////////////////////////////////////////////////////////////////////
// Protected & Private Data
private:
    
    double value;
    Int * owner;
    endpointType endpoint_type;
};

///////////////////////////////////////////////////////////////////////////////
// EndPoint_less
struct EndPoint_less
{
    bool operator()(const EndPoint * p1, const EndPoint * p2) const
    {
        if( (*p1)==(*p2) ){
			if(p1->getType()==p2->getType()) return false; //exact the same...
            if(p1->getType()==EndPoint::LEFT_ENDPOINT) return true;
            if(p2->getType()==EndPoint::LEFT_ENDPOINT) return false;	
        }
        return (*p1)<(*p2);
    }
};

ostream & operator << (ostream & out, const EndPoint& endpoint);

#endif //_ITREE_ENDPOINT_H_
