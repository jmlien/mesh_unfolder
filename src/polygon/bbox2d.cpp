#include "bbox2d.h"
#include "chull.h"
#include "simple_svg_1.0.0.hpp"

namespace masc {
namespace polygon {

//initialize with a given polygon
bbox2d::bbox2d(const c_polygon & poly)
#if DEBUG
: m_chull_ply(c_ply::POUT), m_ply(poly)
#endif
{
  const c_ply& ply=poly.front();
  list<ply_vertex*> hull;
  hull2d(ply.getHead(), ply.getHead()->getPre(), hull);

#if DEBUG
  m_chull_ply.beginPoly();
#endif
  for(ply_vertex* v:hull)
  {
    auto pos=v->getPos();
    this->m_chull.push_back(pos);
#if DEBUG
    m_chull_ply.addVertex(pos[0],pos[1]);
#endif
  }//end for

#if DEBUG
  m_chull_ply.endPoly();
#endif
}

obb bbox2d::build(bbox2d_problem & problem)
{
  mathtool::Vector2d v,n;
  int e[4]; //vertex indices of extreme points
  float a[4]; //angles
  int w; //index (0~3), so that e[w] has the smallest value a[w]
  int hullsize=m_chull.size();

  //init extreme points
  e[0]=0;
  v=(m_chull[1]-m_chull[0]).normalize();
  n.set(-v[1],v[0]);
  const mathtool::Vector2d v0=v;

  float max_v=-FLT_MAX, min_v=FLT_MAX, max_n=-FLT_MAX;
  for(int i=2;i<hullsize;i++)
  {
    auto& pt=m_chull[i];
    double dv=(pt-m_chull[0])*v;
    double dn=(pt-m_chull[0])*n;
    if(dv>max_v){ max_v=dv; e[1]=i;}
    if(dv<min_v){ min_v=dv; e[3]=i;}
    if(dn>max_n){ max_n=dn; e[2]=i;}
  }
  w=findAngles(e,a,v,n);

  //update extreme points
  char svg_filename[256];

  for(int i=0;i<m_chull.size();i++)
  {
    //create a box from v,n,e[4]
    obb box=createOBB(e,v,n);

#if DEBUG
    cout<<"box="<<box<<endl;
    sprintf(svg_filename, "%s%03d.svg", "DEBUG_",i);
    saveSVG(svg_filename,m_chull_ply,box);
    //saveSVG(svg_filename,m_ply.front(),box);
#endif

    //check if this box solve the problem
    if(problem.solved(box)) break;

    //update
    int ne=(e[w]+1)%hullsize;
    mathtool::Vector2d nev=(m_chull[ne]-m_chull[e[w]]).normalize();
    if(w==0 || w==2)
    {
      v=nev;
      n.set(-v[1],v[0]);
    }
    else{
      n=nev;
      v.set(-n[1],n[0]);
    }
    e[w]=ne;

    w=findAngles(e,a,v,n);
  }

  return problem.getSolution(); //done
}

int bbox2d::findAngles
(int e[4], float a[4], const mathtool::Vector2d& v, const mathtool::Vector2d& n)
{
  int size=m_chull.size();
  mathtool::Vector2d u[4];
  for(int i=0;i<4;i++)
    u[i]=(m_chull[(e[i]+1)%size]-m_chull[e[i]]).normalize();

  int w=0;
  a[0]=fabs(v*u[0]);
  a[1]=fabs(n*u[1]); if(a[1]>a[w]){ w=1; } //larger dot product means smaller angle
  a[2]=fabs(v*u[2]); if(a[2]>a[w]){ w=2; }
  a[3]=fabs(n*u[3]); if(a[3]>a[w]){ w=3; }

  // cout<<"e=";
  // for(int i=0;i<4;i++) cout<<e[i]<<", ";
  // cout<<endl;
  //
  // cout<<"a=";
  // for(int i=0;i<4;i++) cout<<a[i]<<", ";
  // cout<<endl;

  return w;
}

obb bbox2d::createOBB(int e[4],const mathtool::Vector2d& v, const mathtool::Vector2d& n)
{
  obb box;
  box.corners[0]=m_chull[e[0]] + v*((m_chull[e[1]]-m_chull[e[0]])*v);
  box.corners[3]=m_chull[e[0]] + v*((m_chull[e[3]]-m_chull[e[0]])*v);
  box.corners[1]=m_chull[e[2]] + v*((m_chull[e[1]]-m_chull[e[2]])*v);
  box.corners[2]=m_chull[e[2]] + v*((m_chull[e[3]]-m_chull[e[2]])*v);

  box.width  = fabs((m_chull[e[3]]-m_chull[e[1]])*v);
  box.height = fabs((m_chull[e[2]]-m_chull[e[0]])*n);

  return box;
}

//the problem of finding the minimum area bounding box
bool min_area_bbox::solved(const obb& box)
{
  float area=box.width*box.height;
  if(area<m_min_area){
    m_min_area=area;
    this->m_solution=box;
  }
  return false; //always return false so search continues;
}

//the problem of finding the minimum boundary bounding box
bool min_perimeter_bbox::solved(const obb& box)
{
  float peri=box.width+box.height;
  if(peri<m_min_peri){
    m_min_peri=peri;
    this->m_solution=box;
  }
  return false; //always return false so search continues;
}

//the problem of finding a bounding box that can fit into another box
bool contained_bbox::solved(const obb& box)
{
  if( (box.width<=m_width && box.height <= m_height) ||
      (box.height<=m_width && box.width <= m_height) )
  {
    this->m_solution=box;
    return true;
  }

  return false;
}

//stream out the box
std::ostream & operator<<(std::ostream& out, const obb& box)
{
  out<<"[w="<<box.width<<", h="<<box.height<<"], ("
     <<box.corners[0]<<") ,("<<box.corners[1]<<"), ("
     <<box.corners[2]<<") ,("<<box.corners[3]<<")";
  return out;
}

//function for saving svg file

void ply2ply(const masc::polygon::c_ply& ply, svg::Polygon& poly)
{
    auto v=ply.getHead();
    do{
      auto & pos=v->getPos();
      poly << svg::Point(pos[0], pos[1]);
      v=v->getNext();
    }
    while(v!=ply.getHead());
    poly.endBoundary();
}

void box2ply(const masc::polygon::obb& box, svg::Polygon& poly)
{
    poly << svg::Point(box.corners[0][0], box.corners[0][1]);
    poly << svg::Point(box.corners[1][0], box.corners[1][1]);
    poly << svg::Point(box.corners[2][0], box.corners[2][1]);
    poly << svg::Point(box.corners[3][0], box.corners[3][1]);
    poly.endBoundary();
}

void saveSVG(string svg_filename, masc::polygon::c_ply& ply, const masc::polygon::obb& box)
{
  //cout<<"??"<<endl;
    //create a svg file
    //char svg_filename[256];
    //sprintf(svg_filename, "%s.svg", img_name.c_str());
    auto R=ply.getRadius();
    auto center=ply.getCenter();

    //cout<<"??"<<endl;
    svg::Dimensions dimensions(R*2.5, R*2.5);
    svg::Document doc(svg_filename, svg::Layout(dimensions, svg::Layout::BottomLeft, 1, svg::Point(-center[0]+R*1.25, -center[1]+R*1.25)));

    //------------------------------------------------------------------
    //draw the external boundary
    ///cout<<"boundary"<<endl;
    svg::Polygon box_bd(svg::Fill(svg::Color::Yellow), svg::Stroke(0.5, svg::Color::Black));
    box2ply(box, box_bd);
    doc << box_bd;

//cout<<"poly"<<endl;
    svg::Polygon poly_bd(svg::Fill(svg::Color::Silver), svg::Stroke(0.5, svg::Color::Black));
    ply2ply(ply, poly_bd);
    doc << poly_bd;

    doc.save();
    cout << "- Saved " << svg_filename << endl;
}

}}//end namespaces
