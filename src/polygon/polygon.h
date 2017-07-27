//------------------------------------------------------------------------------
//  Copyright 2010-2015 by Jyh-Ming Lien and George Mason University
//  See the file "LICENSE" for more information
//------------------------------------------------------------------------------

#pragma once

#ifdef WIN32
#pragma warning(disable : 4786)
#endif

#include "mathtool/Point.h"
#include "mathtool/Vector.h"
using namespace mathtool;

#include <list>
#include <cassert>
#include <vector>
using namespace std;

#include <limits.h>
#include <float.h>


#include "diagonal2.h"

namespace masc {
namespace polygon {

typedef unsigned int uint;

//
//a triangle
//
struct triangle
{
    uint v[3]; // id to the vertices
};


//
// Extra information for Vertex of polygon
//

class c_BPC; //defined in bpc.h, a class for bridge, pocket and concavity

class ply_vertex;
class c_ply;
class c_polygon;
struct ply_vertex_extra
{

    ply_vertex_extra()
    {
		  reInit();
	  }

    void reInit()
    {
    	concavity_bpc=NULL;
      concavity=0;
      flag=0;
      other_v=NULL;
    	accLen=0;
      mesh_vid = UINT_MAX;
      mesh_eid = UINT_MAX;
      on_hull=false;
    }

    //is this vertex a pocket minimum
    bool isPM() const { return (concavity_bpc!=NULL && concavity>0); }

	  ply_vertex * getDihedralPre();  //get previous vertex for estimating the dihedral angle
    ply_vertex * getDihedralNext(); //get next vertex for estimating the dihedral angle

    c_BPC * concavity_bpc; //this bpc defines the concavity of this vertex

    float   concavity;
    uint    flag;
    vector<c_diagonal> diagonals;

    ply_vertex * other_v; //copied vertex in original/simplified polygon

    static uint getFlagID(){ static uint id=1; return id++; }

	  float accLen;		//accumulated length from start vertex... Added by Zhonghua 04/06/2013

	  uint mesh_vid; //the id of mesh vertex that create this ply vertex
	  uint mesh_eid; //the edge id the mesh between this vertex and next vertex
    bool on_hull; //if this vertex is on convex hull
};


//
// Vertex of polygon
//
class ply_vertex
{
public:

    ///////////////////////////////////////////////////////////////////////////
    ply_vertex(){ init(); }
    ply_vertex( const Point2d& p ){ pos=p; init(); }
    virtual ~ply_vertex();
    void setNext(ply_vertex * n){next=n; if(n!=NULL) n->pre=this; }
    void setPre(ply_vertex * n){pre=n; if(n!=NULL) n->next=this; }
    void computeExtraInfo();

    //negate the vertex
    void negate();

    //reverse the order
    void reverse();

    //copy
    void copy(ply_vertex * other);

    ///////////////////////////////////////////////////////////////////////////
    void setPos(const Point2d& p) { pos=p; }
    virtual const Point2d& getPos() const { return pos; }
	  float distanceTo(ply_vertex* other);

    void translate(const Vector2d& v){ pos=pos+v; }

    void rotate(double r);

    virtual ply_vertex * getNext() const { return next; }
    virtual ply_vertex * getPre() const { return pre; }

    const Vector2d& getNormal() const { return normal; }
    bool isReflex() const { return reflex; }

    //get extra information
    uint getVID() const { return vid; }
    void setVID(uint id) {vid=id;}
    ply_vertex_extra& getExtra() { return extra; }
    const ply_vertex_extra& getExtra() const { return extra; }


private:

    //extra info for decomposition
    ply_vertex_extra extra;

    void init(){
        next=pre=NULL;
        reflex=false;
        vid=UINT_MAX;
    }

    //basic info
    Point2d pos;       //position
    ply_vertex * next; //next vertex in the polygon
    ply_vertex * pre;  //previous vertex in the polygon
    Vector2d normal;   //normal, the segment normal from this v to the next.
    bool reflex;
    uint vid;
};

//
// Polygon chain
//
class c_ply{
public:

    enum POLYTYPE { UNKNOWN, PIN, POUT };

    ///////////////////////////////////////////////////////////////////////////
    c_ply(POLYTYPE t){ head=tail=NULL; type=t; radius=-1; area=arclength=-FLT_MAX; canBeIgnored = false;}

    ///////////////////////////////////////////////////////////////////////////
    void copy(const c_ply& ply); //copy from the other ply
    void destroy();

    ///////////////////////////////////////////////////////////////////////////
    // create c_ply
    void beginPoly();
    ply_vertex * addVertex( double x, double y, bool remove_duplicate=false );
    ply_vertex * addVertex( ply_vertex * v );
    void endPoly(bool remove_duplicate=false);

    ///////////////////////////////////////////////////////////////////////////
    void negate();
    void reverse(); //reverse vertex order
    void reverseType(); //revise pin to pout or vice versa

    ///////////////////////////////////////////////////////////////////////////
    void translate(const Vector2d& v);
    void rotate(double radius);
    void scale(float f);

    ///////////////////////////////////////////////////////////////////////////
    // Access
    ply_vertex * getHead() const { return head; }
    POLYTYPE getType() const { return type; }
    void set(POLYTYPE t,ply_vertex * h){
        type=t; head=h;
        if(h!=NULL){ tail=h->getPre(); }
        else{ tail=NULL; }
    }
    int getSize() {
        if(all.empty()) build_all();
        return all.size();
    }

    ply_vertex * operator[](unsigned int id){
        if(all.empty()) build_all();
        return all[id];
    }

    ///////////////////////////////////////////////////////////////////////////
    // additional functions
    const Point2d& getCenter();

    ///////////////////////////////////////////////////////////////////////////
    //compute the Radius of the poly chain
    float getRadius();

    //arc legnth
    float getArcLength();

    //area
    float getArea();

    //check if convex
    bool is_convex() const;

    //delete a vertex
    void delete_vertex(ply_vertex * p);

    ///////////////////////////////////////////////////////////////////////////
    // Operator
    //check if give poly line is the same as this
    bool operator==(const c_ply& other ) const{ return other.head==head; }
    friend istream& operator>>( istream&, c_ply& );
    friend ostream& operator<<( ostream&, c_ply& );

    bool canBeIgnored;

    ///////////////////////////////////////////////////////////////////////////
    bool doInit(); /*return # of vertice in this poly*/

    //build elements in vector<ply_vertex*> all
    void build_all();

protected:



private:

    ply_vertex * head; //the head of vertex list
    ply_vertex * tail; //end of the vertex list

    vector<ply_vertex*> all; //all vertices

    //additional info
    Point2d center;
    float radius;
    float area;
    float arclength;

    //In, out or unknown.
    POLYTYPE type;
};


//a c_plylist is a list of c_ply
class c_plylist : public list<c_ply>
{
    friend ostream& operator<<( ostream&, c_plylist& );
    friend istream& operator>>( istream&, c_plylist& );

public:

    c_plylist()
    {
        box[0]=box[1]=box[2]=box[3]=0;
        is_buildboxandcenter_called=false;
    }

    void negate();
    void translate(const Vector2d& v);
    void rotate(double r);

    //access
    void buildBoxAndCenter();
    double * getBBox() { assert(is_buildboxandcenter_called); return box; }
    const Point2d& getCenter() { assert(is_buildboxandcenter_called); return center; }

protected:

    Point2d center;
    double box[4];

private:

    bool is_buildboxandcenter_called;
};

//
// a c_polygon is a restricted kind of c_plylist
// this defines a simple polygon so that
// the first element much be a POUT c_ply and
// the rest ply lines are a list of holes
//
class c_polygon : public c_plylist
{
public:

    c_polygon() { area=arclength=0; bCnveXEnough = false;}

    void push_back(const c_ply& ply)
    {
        c_plylist::push_back(ply);
        all.clear();
        build_all();
    }

    bool valid(); //check if this is a valid polygon

    //copy from the given polygon
    void copy(const c_polygon& other);

    list<c_polygon> split();

    void reverse(); //reverse the vertex order (not the list order)

    void scale(float factor);

    void normalize();

    //access the vertices of the polygon as an array
    uint getSize()
    {
        if(all.empty()) build_all();
        return all.size();
    }

    //get number of vertices
    uint getSize() const
    {
        assert(all.empty()==false);
        return all.size();
    }


    ply_vertex * operator[](unsigned int id){
        if(all.empty()) build_all();
        return all[id];
    }

    ply_vertex * operator[](unsigned int id) const {
        assert(all.empty()==false);
        return all[id];
    }

    float getArea();

    float getArcLength();

    //destroy
    void destroy();

    bool is_convex() const;

    void build_all();

    bool bCnveXEnough;

private:


    vector<ply_vertex*> all; //all vertices

    float area;
    float arclength;
};

}//end namespace polygon
}//end namespace masc
