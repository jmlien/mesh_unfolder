
//------------------------------------------------------------------------------
//  Copyright 2007-2011 by Jyh-Ming Lien and George Mason University
//  See the file "LICENSE" for more information
//------------------------------------------------------------------------------

#include "polygon.h"
#include "intersection.h"
#include "bpc.h"
#include <vector>


namespace masc {
namespace polygon {

ply_vertex * ply_vertex_extra::getDihedralPre()
{
    assert(other_v);

	//comment by Guilin
    //assert(isPM());
    return other_v->getPre()->getExtra().other_v;
}

ply_vertex * ply_vertex_extra::getDihedralNext()
{
    assert(other_v);

    //comment by Guilin
    //assert(isPM());
    return other_v->getNext()->getExtra().other_v;
}


//
//ply_vertex * ply_vertex_extra::getDihedralPre(){
//    assert(concavity_bpc);
//    ply_vertex * v=concavity_bpc->getConcavity();
//    for(short i=0;i<4;i++) v=v->getPre();
//    return v;
//}
//ply_vertex * ply_vertex_extra::getDihedralNext(){
//    assert(concavity_bpc);
//    ply_vertex * v=concavity_bpc->getConcavity();
//    for(short i=0;i<4;i++) v=v->getNext();
//    return v;
//}


ply_vertex::~ply_vertex()
{
    //doing nothing for now
}

// - compute normal
// - check if the vertex is reflex or not
void ply_vertex::computeExtraInfo()
{
    reflex=true;

    //compute normal direction
    if(next!=NULL){
        Vector2d v=next->pos-pos;
        if( v[0]==0 ){
            if(v[1]>0){ normal[0]=1; normal[1]=0; }
            else{ normal[0]=-1; normal[1]=0; }
        }
        else if( v[0]>0 ){
            normal[1]=-1;
            normal[0]=(v[1]/v[0]);
        }
        else{//v[0]<0
            normal[1]=1;
            normal[0]=-(v[1]/v[0]);
        }

        normal=normal.normalize();

        if(pre!=NULL){
            //compute if left or right turn
            Vector2d u=pos-pre->pos;
            float z=u[0]*v[1]-u[1]*v[0];
            if(z>0) reflex=false;
        }
    }//end if(next!=NULL)
}

void ply_vertex::negate()
{
    normal=-normal;
    pos[0]=-pos[0];
    pos[1]=-pos[1];
}

void ply_vertex::reverse()
{
    swap(next,pre);
    //normal=-normal;
    computeExtraInfo();
}

void ply_vertex::copy(ply_vertex * other)
{
    pos=other->pos;
    normal=other->normal;
    reflex=other->reflex;
    vid=other->vid;
}

void ply_vertex::rotate(double r)
{
    double cos_r=cos(r);
    double sin_r=sin(r);
    //rotate pos
    double x=pos[0]*cos_r-pos[1]*sin_r;
    double y=pos[0]*sin_r+pos[1]*cos_r;
    pos.set(x,y);

    //rotate normal
    x=normal[0]*cos_r-normal[1]*sin_r;
    y=normal[0]*sin_r+normal[1]*cos_r;
    normal.set(x,y);
}

float ply_vertex::distanceTo(ply_vertex* other)
{
	if(other == NULL) assert(false);

	return (this->pos - other->pos).norm();
}

///////////////////////////////////////////////////////////////////////////////

//copy from the given ply
void c_ply::copy(const c_ply& other)
{
    destroy();//detroy myself first

    ply_vertex* ptr=other.head;
    beginPoly();
    do{
        ply_vertex * v=new ply_vertex();
        assert(v); //check for memory
        v->copy(ptr);
        addVertex(v);
        ptr=ptr->getNext();
    }while( ptr!=other.head );

    //endPoly();
    //finish up
    tail->setNext(head);

    //copy extra info
    area=other.area;
    center=other.center;
    radius=other.radius;
    type=other.type;
}

// clean up the space allocated
void c_ply::destroy()
{
    if( head==NULL ) return;
    ply_vertex* ptr=head;
    do{
        ply_vertex * n=ptr->getNext();
        delete ptr;
        ptr=n;
    }
    while( ptr!=head && ptr!=NULL);
    head=tail=NULL;

    all.clear();
}

// Create a empty polygon
void c_ply::beginPoly()
{
    head=tail=NULL;
    all.clear();
}

// Add a vertex to the polygonal chian
ply_vertex * c_ply::addVertex( double x, double y, bool remove_duplicate )
{
    Point2d pt(x,y);

    if(tail!=NULL){
        if(tail->getPos()==pt && remove_duplicate){
            cout<<"duplicated vertex"<<endl;
            return NULL; //don't add
        }
    }

    ply_vertex * v=new ply_vertex(pt);
    if( tail!=NULL ){
        tail->setNext(v);
    }
    tail=v;
    if( head==NULL ) head=tail;
    v->setVID(all.size()); //id of the vertex in this ply
    all.push_back(v);

    return v;
}

// Add a vertex to the polygonal chian
ply_vertex * c_ply::addVertex( ply_vertex * v )
{
    if( tail!=NULL ){
        tail->setNext(v);
    }
    tail=v;
    if( head==NULL ) head=tail;
    v->setVID(all.size()); //id of the vertex in this ply
    all.push_back(v);

    return v;
}

// finish building the polygon
void c_ply::endPoly(bool remove_duplicate)
{
    if(head!=NULL && tail!=NULL){
        if(remove_duplicate){
            if(head->getPos()==tail->getPos()){ //remove tail..
                delete tail;
                all.pop_back();
                tail=all.back();
            }
        }//
    }

    tail->setNext(head);
    doInit();
}

// initialize property of the this polychain
// Compute normals and find reflective vertices
bool c_ply::doInit()
{
	bool b = true;
    //compute area
    getArea();
    if(this->area<0 && type==POUT){
       //cerr<<"! Warning: polygon type is POUT but has negative area. Reverse the vertex ordering."<<endl;
       reverse();

       b = false;
    }
    else if(this->area>0 && type==PIN){
       //cerr<<"! Warning: polygon type is PIN but has positive area. Reverse the vertex ordering."<<endl;
       reverse();

       b = false;
    }

    //compute normals
    ply_vertex* ptr=head;
    do{
        ptr->computeExtraInfo();
        ptr=ptr->getNext();
    }while( ptr!=head );

    return b;
}

const Point2d& c_ply::getCenter()
{
    if(radius<0){
        center.set(0,0);
        ply_vertex * ptr=head;
        const Point2d& first=ptr->getPos();
        uint size=0;
        do{
            size++;
            Vector2d v=ptr->getPos()-first;
            center[0]+=v[0];
            center[1]+=v[1];
            ptr=ptr->getNext();
        }while(ptr!=head); //end while
        center[0]=(center[0]/size)+first[0];
        center[1]=(center[1]/size)+first[1];

        radius=0;
    }

    return center;
}


///////////////////////////////////////////////////////////////////////////
void c_ply::negate()
{
    ply_vertex * ptr=head;
    do{
        ptr->negate();
        ptr=ptr->getNext();
    }while(ptr!=head); //end while
}

//reverse the order of the vertices
void c_ply::reverse()
{
    ply_vertex * ptr=head;
    do{
        ptr->reverse();
        ptr=ptr->getNext();
    }
    while(ptr!=head); //end while

    this->area=-this->area;
    all.clear();
}

void c_ply::reverseType()
{
    reverse();
    type=(type==PIN)?POUT:PIN;
}

///////////////////////////////////////////////////////////////////////////
void c_ply::translate(const Vector2d& v)
{
    ply_vertex * ptr=head;
    do{
        ptr->translate(v);
        ptr=ptr->getNext();
    }while(ptr!=head); //end while
}

void c_ply::rotate(double radius)
{
    ply_vertex * ptr=head;
    do{
        ptr->rotate(radius);
        ptr=ptr->getNext();
    }while(ptr!=head); //end while
}


///////////////////////////////////////////////////////////////////////////
//compute the Radius of the poly chain
float c_ply::getRadius()
{
    if(radius<0) getCenter();

    if(radius==0){
        ply_vertex * ptr=head;
        do{
            float d=(center-ptr->getPos()).normsqr();
            if(d>radius) radius=d;
            ptr=ptr->getNext();
        }while(ptr!=head); //end while
        radius=sqrt(radius);
    }

    return radius;
}

float c_ply::getArcLength()
{
	if (arclength == -FLT_MAX){
		arclength = 0;
		ply_vertex * ptr = head;
		do{
			ply_vertex * next = ptr->getNext();
			const Point2d& p1 = ptr->getPos();
			const Point2d& p2 = next->getPos();
			arclength += (p1 - p2).norm();
			ptr = next;
		} while (ptr != head); //end while
	}

	return arclength;
}

float c_ply::getArea()
{
    if(area==-FLT_MAX){
        area=0;
        getCenter();
        ply_vertex * ptr=head;
        do{
            ply_vertex * next=ptr->getNext();
            const Point2d& p1=ptr->getPos();
            const Point2d& p2=next->getPos();
            area+=Area(p1.get(),p2.get(),center.get());
            ptr=next;
        }while(ptr!=head); //end while
    }

    return area;
}


//check if convex
bool c_ply::is_convex() const
{
    ply_vertex * ptr=head;
    do{
        if(ptr->isReflex()) return false;
        ptr=ptr->getNext();
    }while(ptr!=head); //end while

    return true;
}

void c_ply::delete_vertex(ply_vertex * v)
{
    ply_vertex *pre=v->getPre();
    ply_vertex *next=v->getNext();

    pre->setNext(next);
    next->setPre(pre);
    pre->computeExtraInfo(); //recompute info

    if(head==v){
        head=next;
        tail=pre;
    }
    delete v;

    all.clear(); //not valid anymore
}

void c_ply::build_all()
{
    uint vid=0;
    all.clear();
    ply_vertex * ptr=head;
    do{
        ptr->setVID(vid++);
        all.push_back(ptr);
        ptr=ptr->getNext();
    }while(ptr!=head); //end while
}

void c_ply::scale(float f)
{
    ply_vertex * ptr=head;
    do{
        const Point2d& p = ptr->getPos();
        Point2d np(p[0]*f,p[1]*f);
        ptr->setPos(np);
        ptr=ptr->getNext();
    }while(ptr!=head); //end while

    area=-FLT_MAX;
    radius=-FLT_MAX;
}

//
// Compute the center and the box of a list of plys
//

void c_plylist::buildBoxAndCenter()
{
    //typedef list<c_ply>::iterator IT;
    box[0]=box[2]=FLT_MAX;
    box[1]=box[3]=-FLT_MAX;
    for(iterator i=begin();i!=end();i++){

        ply_vertex * ptr=i->getHead();
        do{
            const Point2d& p=ptr->getPos();
            if(p[0]<box[0]) box[0]=p[0];
            if(p[0]>box[1]) box[1]=p[0];
            if(p[1]<box[2]) box[2]=p[1];
            if(p[1]>box[3]) box[3]=p[1];
            ptr=ptr->getNext();
        }
        while(ptr!=i->getHead()); //end while
    }

    center[0]=(box[0]+box[1])/2;
    center[1]=(box[2]+box[3])/2;

    is_buildboxandcenter_called=true;
}

//
// Compute the center and the box of a list of plys
//
void c_plylist::translate(const Vector2d& v)
{
    for(iterator i=begin();i!=end();i++) i->translate(v);
}

//
// rotate all polychains
//
void c_plylist::rotate(double r)
{
    for(iterator i=begin();i!=end();i++) i->rotate(r);
}

void c_plylist::negate()
{
    for(iterator i=begin();i!=end();i++) i->negate();
}

void c_polygon::reverse()
{
    for(iterator i=begin();i!=end();i++) i->reverse();
    all.clear();
}

bool c_polygon::valid() //check if this is a valid polygon
{
    typedef list<c_ply>::iterator IT;
    if(empty()) return false;
    if(front().getType()!=c_ply::POUT) return false;
    for(iterator i=++begin();i!=end();i++) if(i->getType()!=c_ply::PIN) return false;

    return true;
}

//copy from the given polygon
void c_polygon::copy(const c_polygon& other)
{
    clear();
    all.clear();

    for(const_iterator i=other.begin();i!=other.end();i++){
        c_ply p(c_ply::UNKNOWN);
        p.copy(*i);
        push_back(p);
    }
}

list<c_polygon> c_polygon::split()
{
    list<c_polygon> sub;
    for(iterator i=begin();i!=end();i++){
        if(i->getType()==c_ply::POUT){
            c_polygon tmp;
            tmp.push_back(*i);
            sub.push_back(tmp);
        }
        else{
            if(sub.empty()){
                cerr<<"! Error: Invalid polygon type: Holes are defined without external boundary."<<endl;
                continue; //can't do anything about this...
            }
            sub.back().push_back(*i);
        }
    }
    return sub;
}

void c_polygon::scale(float factor)
{
    for(iterator i=begin();i!=end();i++)
        i->scale(factor);
}

void c_polygon::normalize()
{
    float r=front().getRadius();
    scale(1.0/r);
}


void c_polygon::destroy()
{
    for(iterator i=begin();i!=end();i++){
        i->destroy();
    }
    clear(); //remove all ply from this list
    all.clear();
}

bool c_polygon::is_convex() const
{
    if(size()>1) return false; //contains hole
    return front().is_convex();
}

void c_polygon::build_all()
{
    uint vid=0;
    for(iterator i=begin();i!=end();i++){
        uint vsize=i->getSize();
        for(uint j=0;j<vsize;j++){
            (*i)[j]->setVID(vid++);       //id of vertex in the polygons
            all.push_back((*i)[j]);
        }
    }//end for i
}

float c_polygon::getArea()
{
    if(area==0){
        for(iterator i=begin();i!=end();i++){
            area+=i->getArea();
        }//end for i
    }

    return area;
}


float c_polygon::getArcLength()
{
	if (arclength == 0){
		for (iterator i = begin(); i != end(); i++){
			arclength += i->getArcLength();
		}//end for i
	}

	return arclength;
}


istream& operator>>( istream& is, c_ply& poly)
{
    int vsize; string str_type;
    is>>vsize>>str_type;

    if( str_type.find("out")!=string::npos )
        poly.type=c_ply::POUT;
    else poly.type=c_ply::PIN;

    poly.beginPoly();
    //read in all the vertices
    int iv;
    vector< pair<double,double> > pts; pts.reserve(vsize);
    for( iv=0;iv<vsize;iv++ ){
        double x,y;
        is>>x>>y;
        pts.push_back(pair<double,double>(x,y));
        //double d=x*x+y*y;
    }
    int id;
    for( iv=0;iv<vsize;iv++ ){
        is>>id; id=id-1;
        poly.addVertex(pts[id].first,pts[id].second,true);
    }

    poly.endPoly();
    return is;
}

istream& operator>>( istream& is, c_plylist& p)
{
    //remove header commnets
    do{
        char tmp[1024];
        char c=is.peek();
        if(isspace(c)) is.get(c); //eat it
        else if(c=='#') {
            is.getline(tmp,1024);
        }
        else break;
    }while(true);

    //start reading
    uint size;
    is>>size;
    uint vid=0;
    for(uint i=0;i<size;i++){
        c_ply poly(c_ply::UNKNOWN);
        is>>poly;
        p.push_back(poly);
        uint vsize=poly.getSize();
        for(uint j=0;j<vsize;j++){
            poly[j]->setVID(vid++);       //id of vertex in the polygons
            //p.all.push_back(poly[j]);
        }
    }
    return is;
}

ostream& operator<<( ostream& os, c_ply& p)
{
    os<<p.getSize()<<" "<<((p.type==c_ply::PIN)?"in":"out")<<"\n";
    ply_vertex * ptr=p.head;
    do{
        os<<ptr->getPos()[0]<<" "<<ptr->getPos()[1]<<"\n";
        ptr=ptr->getNext();
    }while(ptr!=p.head);

    for(int i=0;i<p.getSize();i++) os<<i+1<<" ";
    os<<"\n";
    return os;
}

ostream& operator<<( ostream& out, c_plylist& p)
{
    out<<p.size()<<"\n";
    typedef c_plylist::iterator PIT;
    for(PIT i=p.begin();i!=p.end();i++) out<<*i;
    return out;
}


}//end namespace polygon
}//end namespace masc
