/*
 * SVGWriter.h
 *
 *  Created on: Nov 21, 2016
 *      Author: zxi
 */

#ifndef SRC_UTIL_SVGWRITER_H_
#define SRC_UTIL_SVGWRITER_H_

#include <vector>
#include <list>
#include <string>
#include <iostream>
#include <map>
#include <memory>
using namespace std;


#include "model.h"
#include "config.h"
#include "mathtool/Point.h"
using namespace mathtool;

#include "util/TextureRenderer2D.h"
using namespace masc::util;

#include "intersection.h"

// forward declaration
struct model;

namespace masc {
namespace unfolding {
namespace util {

enum class ExportSVGType {

  CUT = 1,                       // cut, boundary in black, creases in dotted gray lines.
  BASIC = 2,                     // show edges in colors
  TREE_ONLY = 4,                 // only show the tree, private
  FLAT_EDGES_ONLY = 8,           //
  NUMBERS = 16,                  // show number of cut edges
  TEXTURE = 32,                  // show texture
  NORMAL = BASIC | NUMBERS | TEXTURE,      // basic with numbers
  NET = BASIC | FLAT_EDGES_ONLY, // show all edges
  TREE = NET | TREE_ONLY,        // net mode with tree
  ALL_EDGES = NORMAL | FLAT_EDGES_ONLY
//TODO(zxi) add mode
};

inline int operator&(ExportSVGType a, ExportSVGType b) {
  return static_cast<int>(a) & static_cast<int>(b);
}

inline ExportSVGType operator|(ExportSVGType a, ExportSVGType b) {
  return static_cast<ExportSVGType>(static_cast<int>(a) & static_cast<int>(b));
}

/*
 Usage:
 SVGWriter writer(unfolder->getNet(), unfolder->getConfig());
 writer.Init();
 writer.Save("tmp.svg", ExportSVGType::CUT);
 */

class SVGWriter {
public:
    // Create an SVG writer, does not own the model.
    // the model must be flattened onto XZ plane.
    SVGWriter(const model* model, const Config& config) :
      inited_(false), model_(model), config_(config)
    {
        width_pixel_ = height_pixel_ = 0;
        font_size_ = 1.0;
    }

    // Initialize the styles, compute the bounding box, boundary, etc.
    void Init();

    // Get Border Cuts edge list
    list< uint >* GetBorderCuts() {
        this->AnalyzeNet();
        return &this->border_cut_eids_;
    }

    ~SVGWriter() { }

    void Save(const string& output_path, ExportSVGType type);

    // Find the boundary for the given model.
    void FindBoundaryPolygon(vector<int>* boundary_vertices);

    // Get the coordinate for a given vertex in the model.
    Vector3d GetSVGCoord(uint vid) const;

private:

  struct Tab
  {
      bool intersect(const Tab& other) const
      {
        double pp[2] = { 0.0, 0.0 };
        auto it=shape_.begin();
        auto ne=it; ne++;
        for(;ne!=shape_.end();it=ne,ne++)
        {
          const double a1[2] = {(*it)[0], (*it)[2]};
          const double b1[2] = {(*ne)[0], (*ne)[2]};

          auto it2=other.shape_.begin();
          auto ne2=it2; ne2++;
          for(;ne2!=other.shape_.end();it2=ne2,ne2++)
          {
            const double a2[2] = {(*it2)[0], (*it2)[2]};
            const double b2[2] = {(*ne2)[0], (*ne2)[2]};
            if( SegSegInt<double>(a1, b1, a2, b2, pp)=='1' ) { return true;}
          }//end for ne2
        }//end for ne
        return false;
      }

      bool intersect(const Vector3d & u, const Vector3d & v) const
      {

        //regular checks
        double pp[2] = { 0.0, 0.0 };
        const double v2[2] = { v[0], v[2] };
        const double u2[2] = { u[0], u[2] };

        auto it=shape_.begin();
        auto ne=it; ne++;
        for(;ne!=shape_.end();it=ne,ne++)
        {
          double a1[2] = {(*it)[0], (*it)[2]};
          double b1[2] = {(*ne)[0], (*ne)[2]};

          //check if ac intersect uv
          char r = SegSegInt<double>(a1, b1, v2, u2, pp);

          //avoid intersection at end points
          Vector3d tmp(pp[0],0,pp[1]);
          if( (tmp-v).normsqr()<1e-10 ) continue;
          if( (tmp-u).normsqr()<1e-10 ) continue;

          if(r=='1') return true;

        }//end for ne

        return false;
      }

      //Vector3d a, b, c; //ab is on the cut edge
      list<Vector3d> shape_;
      uint eid;         //eid of the edge that this tab is appended to
  };

  void WriteHeader(ostream& out);
  void WriteStyle(ostream& out);
  void WriteTexture(ostream& out);

  void WriteBoundaryPolygon(ostream& out);
  void WriteBoundaryPolygon2(ostream& out);

  void WriteCreases(ostream& out);
  void WriteCreasesInColor(ostream& out);
  void WriteFlatCreases(ostream& out);
  void WriteValleyCreasesOnly(ostream& out);

  void AnalyzeNet();                       //this may change the net boundary
  void WriteCutEdgeID(ostream& out);       //print cut edge ids
  void WriteCutEdgeHints(ostream& out);    //extra hints for connecting cut edges
  void WriteZipHints(ostream& out);        //zip

  void WriteTabs(ostream& out);

  void WriteTree(ostream& out);
  void WriteLabels(ostream & out);
  void WriteFooter(ostream& out);

  // Get the coordinate for a given vertex in the model.
  Vector3d GetSVGCoord(const Point3d& v) const;
  Vector3d GetSVGCoord(const Vector3d& v) const;

  //tab related Methods
  void BuildTabs();
  Tab BuildTab(Vector3d& a, Vector3d& b, Vector3d& dir, double len, uint eid );

  //check if the given tab is valid (ie without intersection with the cut polygon and other tabs)
  bool IsValid(const Tab& tab) const;

  //greedy approach to sort line segments so the total travel distance is minimized
  inline void sort_line_segments(list< pair<Vector3d,Vector3d> >& segs)
  {
    if(segs.empty()) return;

    list< pair<Vector3d,Vector3d> > segs_sorted;
    segs_sorted.push_back(segs.front());
    segs.pop_front();

    while(segs.empty()==false)
    {
      //find a closest seg to the last seg in segs_sorted
      auto& last_seg = segs_sorted.back();
      float min_d=FLT_MAX;
      auto best_seg=segs.end();

      for(auto i=segs.begin();i!=segs.end();i++)
      {
         float d=(i->first-last_seg.second).normsqr();
         if(d<min_d)
         {
           min_d=d;
           best_seg=i;
         }
      }//end for

      assert(best_seg!=segs.end());
      segs_sorted.push_back(*best_seg);
      segs.erase(best_seg);

    }//end while

    segs.swap(segs_sorted);
  }

  inline uint otherv(const triangle& tri, const  edge& e)
  {
    for(uint d=0;d<3;d++)
    {
      if(e.vid[0]!=tri.v[d] && e.vid[1]!=tri.v[d]) return tri.v[d];
    }
    cerr<<"! Error: otherv in SVGWriter failed"<<endl;
    assert(false);
    return UINT_MAX;
  }

  inline uint getrootid(const model* model_, uint vid)
  {
    while( model_->vertices[vid].parent_id != UINT_MAX)
      vid=model_->vertices[vid].parent_id;
    return vid;
  }

  inline uint geteid(uint v1, uint v2) const
  {
    for(uint eid : model_->vertices[v1].m_e)
    {
      const edge& e = model_->edges[eid];
      if( (e.vid[0]==v1 && e.vid[1]==v2) || (e.vid[0]==v2 && e.vid[1]==v1) )
        return eid;
    }
    return UINT_MAX;
  }

  bool inited_;

  // Do not own the objective, must be alive during the life-cycle of the SVGWriter object.
  const model* model_;

  // Do not own the objective.
  const Config& config_;

  vector<int> boundary_vertices_;

  // class_name -> style.
  map<string, string> styles_;

  // Dimension of the net.
  Vector3d dim_;

  // Min point in the file
  Point3d net_min_;

  int width_pixel_;
  int height_pixel_;
  float font_size_;

  Vector3d PADDING;

  std::string texture_path_;
  std::unique_ptr<TextureRenderer2D> texture_render_;

  unordered_map<uint, Tab> tabs_; //the first uint is the eid
  list<uint> border_cut_eids_;
  list< list<uint> > zip_line_vids_;
  list<uint> zip_line_eids_;
};

} /* namespace util */
} /* namespace unfolding */
} /* namespace masc */

#endif /* SRC_UTIL_SVGWRITER_H_ */
