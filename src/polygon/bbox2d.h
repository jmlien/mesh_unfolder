#pragma once

#include "Point.h"
#include "polygon.h"
#include <ostream>

#define DEBUG 0 //to enable debugging output, change to 1 

namespace masc {
namespace polygon {

  //oritented bounding box (obb)
  struct obb
  {
      obb(){ width=height=FLT_MAX; }
      mathtool::Point2d corners[4];
      float width, height;
  };

  class bbox2d_problem
  {
    public:
      //return true if the problem is solved by the given box
      virtual bool solved(const obb& box)=0;
      const obb & getSolution() const { return m_solution; }

    protected:
      obb m_solution;
  };

  class bbox2d
  {
  public:

      //initialize with a given polygon
      bbox2d(const c_polygon & poly);

      //build an obb of the provided polygon that solved the given problem.
      //return: obb that solves the given problem.
      obb build(bbox2d_problem & problem);

  private:

      //return the index of the smallest angle
      int findAngles(int e[4], float a[4], const mathtool::Vector2d& v, const mathtool::Vector2d& n);

      //create create OBB
      obb createOBB(int e[4],const mathtool::Vector2d& v, const mathtool::Vector2d& n);

      //data
#if DEBUG
      c_ply m_chull_ply;
      c_polygon m_ply;
#endif
      vector<mathtool::Point2d> m_chull; //convex hull of the input poly
  };

  //we now define problems of finding various boudning boxes

  //the problem of finding the minimum area bounding box
  class min_area_bbox : public bbox2d_problem
  {
    public:
      min_area_bbox(){ m_min_area=FLT_MAX; }
      bool solved(const obb& box);
    private:
      float m_min_area;
  };

  //the problem of finding the minimum boundary bounding box
  class min_perimeter_bbox : public bbox2d_problem
  {
    public:
      min_perimeter_bbox(){ m_min_peri=FLT_MAX; }
      bool solved(const obb& box);
    private:
      float m_min_peri;
  };

  //the problem of finding a bounding box that can fit into another box
  class contained_bbox : public bbox2d_problem
  {
    public:
      contained_bbox(float width, float height){m_width=width; m_height=height;}
      bool solved(const obb& box);
    private:
      float m_width, m_height;
  };

  //stream out the box
  std::ostream & operator<<(std::ostream& out, const obb& box);

  //save to svg file
  void saveSVG(string svg_filename, c_ply& ply, const obb& box);

}//end namespace polygon
}//end namespace masc
