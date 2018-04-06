/*
 * SVGWriter.cpp
 *
 *  Created on: Nov 21, 2016
 *      Author: zxi
 */

#include "util/SVGWriter.h"
#include "polygon/net_analyzer.h"
#include "polygon/polygon.h"
#include "mathtool/Box.h"

#include <cstdlib>
#include <sstream>

namespace {
const int SVG_PADDING_PIXEL = 8;
}

namespace masc {
namespace unfolding {
namespace util {

//get a Bezier curve from 3 control points a, b,c
//implemented in this file
void Bezier(const Vector3d & a, const Vector3d & b, const Vector3d & c,
    list<Vector3d>& curve, int size);

void SVGWriter::Init() {

  if (inited_)
    return;
  inited_ = true;

  PADDING = Vector3d(SVG_PADDING_PIXEL, SVG_PADDING_PIXEL, SVG_PADDING_PIXEL);

  // Find boundary vertices
  if (boundary_vertices_.empty()) {
    this->FindBoundaryPolygon(&boundary_vertices_);
  }

  // Compute bounding box
  Box3d box;
  for (const vertex& v : model_->vertices) {
    box.addPoint(v.p);
  }

  //expand the box to accomondate tabs
  if (this->config_.add_tabs) {
    Vector3d dim = box.getDim() * 0.55f; //scale must be >0.5
    Vector3d mid = (box.getMin() + box.getMax()) * 0.5;
    box.addPoint(mid + dim);
    box.addPoint(mid - dim);
  }

  dim_ = box.getDim();
  net_min_ = box.getMin();

  const auto width = dim_[0] * config_.scale;
  const auto height = dim_[2] * config_.scale;
  const auto aspect = width / height;

  width_pixel_ = width + 2 * SVG_PADDING_PIXEL;
  height_pixel_ = width_pixel_ / aspect + 2 * SVG_PADDING_PIXEL;

  const auto dashed_length = width_pixel_ * 0.005;

  const auto stroke_width = 1.0 * max(width_pixel_, height_pixel_) / 800;

  font_size_ = min(width_pixel_, height_pixel_) * 0.02
      * config_.label_font_scale;

  const string dashed_line = "stroke-dasharray:" + std::to_string(dashed_length)
      + ", " + std::to_string(dashed_length) + ";";

  const string boundary_style = "stroke:rgb(0,0,0);stroke-width:"
      + std::to_string(stroke_width * 2) + ";fill:none";
  const string normal_style = "stroke:rgb(128,128,128);stroke-width:"
      + std::to_string(stroke_width);
  const string mountain_style = "stroke:rgb(255,0,0);stroke-width:"
      + std::to_string(stroke_width) + ";"; // + dashed_line;

  const string valley_style = "stroke:rgb(0,0,255);stroke-width:"
      + std::to_string(stroke_width) + ";";

  const string hint_style = "stroke:rgb(240,200,0);stroke-width:"
      + std::to_string(stroke_width * 2) + ";";

  const string zip_style = "stroke:rgb(0,240,200);stroke-width:"
      + std::to_string(stroke_width / 2) + ";";

  const string tree_style = "stroke:rgb(0,128,0);stroke-width:"
      + std::to_string(stroke_width) + ";";

  const string crease_style = "stroke:rgb(128,128,128);stroke-width:"
      + std::to_string(stroke_width) + ";" + dashed_line;

  const string cut_style = "stroke:rgb(0,0,0);stroke-width:"
      + std::to_string(stroke_width) + ";" + dashed_line;

  const string tab_style = "stroke:rgb(0,0,255);stroke-width:"
      + std::to_string(stroke_width / 2) + ";fill:none";

  styles_["m"] = mountain_style;
  styles_["v"] = valley_style;
  styles_["b"] = boundary_style;
  styles_["n"] = normal_style;
  styles_["t"] = tree_style;
  styles_["c"] = crease_style;
  styles_["h"] = hint_style;
  styles_["z"] = zip_style;
  styles_["tab"] = tab_style;

  const auto scale_factor = config_.scale;

  // remember texture path
  texture_path_=model_->texture_path;

}

void SVGWriter::Save(const string& output_path, ExportSVGType type) {

  if (!inited_) {
    this->Init();
  }
  //
  ofstream out(output_path);

  if (!out.good()) {
    cerr << "!Error! Failed to open file = " << output_path << endl;
    return;
  }

  this->WriteHeader(out);
  this->WriteStyle(out);

  // Only writes texture in normal mode.
  if (!this->texture_path_.empty() && (type & ExportSVGType::TEXTURE) )
  {
    this->texture_render_.reset(new TextureRenderer2D(this->texture_path_));
    assert(this->texture_render_!=nullptr);
    this->texture_render_->SetCanvasSize(this->width_pixel_, this->height_pixel_);
    this->WriteTexture(out);
  }

  AnalyzeNet();
  if (config_.add_tabs)
    BuildTabs();

  //don't write anything else if the texture is the only thing we need
  if(type != ExportSVGType::TEXTURE)
  {
    this->WriteBoundaryPolygon(out);

    if (type == ExportSVGType::CUT) {
      this->WriteCreases(out);
      if (this->config_.svg_edge_hints) {
        //zip line
        this->WriteZipHints(out);
        //mark valley fold as well
        this->WriteValleyCreasesOnly(out);
        //mark cut edges
        this->WriteCutEdgeHints(out);
        //tabs
        if (config_.add_tabs)
          WriteTabs(out);
      }
    }

    if (type & ExportSVGType::FLAT_EDGES_ONLY) {
      this->WriteFlatCreases(out);
    }

    if (type & ExportSVGType::TREE_ONLY) {
      this->WriteTree(out);
    }

    if (type & ExportSVGType::NORMAL) {
      this->WriteCreasesInColor(out);
    }

    if (type & ExportSVGType::NUMBERS) {
      this->WriteLabels(out);
    }
  }//

  this->WriteFooter(out);

  out.close();

  cerr << "- Wrote SVG to " << output_path << endl;
}

///////////////////////////////////////////////////////////////////////////
// Private Methods, clients can stop reading.
///////////////////////////////////////////////////////////////////////////
void SVGWriter::WriteHeader(ostream& out) {
  out << "<?xml version=\"1.0\"?>" << endl;
  out << "<svg width=\"" << width_pixel_ << "\" height=\"" << height_pixel_
      << "\"" << " version=\"1.1\" xmlns=\"http://www.w3.org/2000/svg\" "
      << "xmlns:xlink=\"http://www.w3.org/1999/xlink\">" << endl;
}

void SVGWriter::WriteStyle(ostream& out) {
  out << "  <defs>" << endl;
  out << "  <style type=\"text/css\">" << endl;
  out << "    <![CDATA[" << endl;

  for (const auto& kv : styles_) {
    out << "      ." << kv.first << " {" << kv.second << "}" << endl;
  }

  out << "    ]]>" << endl;
  out << "  </style>" << endl;
  out << "  </defs>" << endl;
}

void SVGWriter::WriteTexture(ostream& out)
{
  cerr << "- Writing Texture" << endl;

  // For each face
  for (int fid = 0; fid < model_->t_size; ++fid)
  {
    const auto& face = this->model_->tris[fid];

    const Vector3d& p0_3d = GetSVGCoord(face.v[0]);
    const Vector3d& p1_3d = GetSVGCoord(face.v[1]);
    const Vector3d& p2_3d = GetSVGCoord(face.v[2]);

    // Coordinates on SVG canvas
    const Vector2d p0(p0_3d[0], p0_3d[2]);
    const Vector2d p1(p1_3d[0], p1_3d[2]);
    const Vector2d p2(p2_3d[0], p2_3d[2]);

    // UV coordinates
    const Vector2d& uv0 = model_->texture_pts[face.vt[0]];
    const Vector2d& uv1 = model_->texture_pts[face.vt[1]];
    const Vector2d& uv2 = model_->texture_pts[face.vt[2]];

    this->texture_render_->Render(p0, p1, p2, uv0, uv1, uv2);

  }

  std::string img_url;
  this->texture_render_->ExportToPngBase64(&img_url);
  out << "<image xlink:href=\"" << img_url << "\" width=\"" << width_pixel_
      << "px\" height=\"" << height_pixel_ << "px\" />\n";
}

void SVGWriter::WriteBoundaryPolygon(ostream& out) {

  // draw creases in dotted grey as a single path
  out << "  <path fill=\"none\" class=\"b\" d=\"";
  auto vsize = boundary_vertices_.size();
  for (auto i = 0; i < vsize; ++i) {

    Vector3d v = this->GetSVGCoord(boundary_vertices_[i]);

    if (i == 0) {
      out << " M ";
    } else {
      out << " L ";
    }

    out << v[0] << " " << v[2];

    //output tab for the current edge, if there is any
    if (i != 0 && tabs_.empty() == false) {
      int j = (i + 1) % vsize;
      uint eid = this->geteid(boundary_vertices_[i], boundary_vertices_[j]);
      //find if there is a tab created for this edge
      //if( find(border_cut_eids_.begin(), border_cut_eids_.end(), eid)!=border_cut_eids_.end())
      //{
      if (tabs_.find(eid) != tabs_.end()) {
        Tab & tab = tabs_[eid];
        const auto& start = tab.shape_.front();
        Vector3d u = this->GetSVGCoord(boundary_vertices_[j]);
        if ((start - v).normsqr() < 1e-10) {
          for (auto it = ++tab.shape_.begin(); it != tab.shape_.end(); it++)
            out << " L " << (*it)[0] << " " << (*it)[2];
        } else if ((start - u).normsqr() < 1e-10) {
          for (auto it = ++tab.shape_.rbegin(); it != tab.shape_.rend(); it++)
            out << " L " << (*it)[0] << " " << (*it)[2];
        } else {
          cerr
              << "! Error: SVGWriter::WriteBoundaryPolygon: When creating a tab. Ignore this tab."
              << endl;
        }
      }
      //}
      //Vector3d v = this->GetSVGCoord(boundary_vertices_[i]);
    }
  }
  out << "\" />" << endl;
}

void SVGWriter::WriteBoundaryPolygon2(ostream& out) {

  // draw creases in dotted grey as a single path
  out << "  <path class=\"b\" d=\"";

  for (int i = 0; i < boundary_vertices_.size() - 1; ++i) {

    Vector3d v1 = this->GetSVGCoord(boundary_vertices_[i]);
    Vector3d v2 = this->GetSVGCoord(boundary_vertices_[i + 1]);

    out << " M " << v1[0] << " " << v1[2];
    out << " L " << v2[0] << " " << v2[2];
  }
  out << "\" />" << endl;
}

void SVGWriter::WriteCreases(ostream& out) {

  list<pair<Vector3d, Vector3d> > creases;

  // draw creases in dotted grey as a single path
  for (int i = 0; i < model_->e_size; ++i) {
    const edge& e = model_->edges[i];
    if (e.type == 'b' || fabs(e.folding_angle) < 1e-3)
      continue;

    Vector3d p1 = GetSVGCoord(e.vid[0]);
    Vector3d p2 = GetSVGCoord(e.vid[1]);
    creases.push_back(make_pair(p1, p2));
  }

  //tab crease
  for (auto& tab : tabs_) {
    creases.push_back(
        make_pair(tab.second.shape_.front(), tab.second.shape_.back()));
  }

  sort_line_segments(creases);

  out << "  <path class=\"c\" d=\"";
  for (auto& crease : creases) {
    const auto& p1 = crease.first;
    const auto& p2 = crease.second;
    out << " M " << p1[0] << " " << p1[2];
    out << " L " << p2[0] << " " << p2[2];
  }
  out << "\" />" << endl;
}

void SVGWriter::WriteValleyCreasesOnly(ostream& out) {
  list<pair<Vector3d, Vector3d> > creases;

  // draw creases in dotted grey as a single path
  for (int i = 0; i < model_->e_size; ++i) {
    const edge& e = model_->edges[i];
    if (e.folding_angle > -1e-3)
      continue;

    const Vector3d p1 = GetSVGCoord(e.vid[0]);
    const Vector3d p2 = GetSVGCoord(e.vid[1]);
    const auto vec = (p2 - p1) / 4;
    creases.push_back(make_pair(p1 + vec, p2 - vec));
  }

  if (creases.empty())
    return;

  sort_line_segments(creases);

  out << "  <path class=\"h\" d=\"";
  for (auto& crease : creases) {
    const auto& p1 = crease.first;
    const auto& p2 = crease.second;
    out << " M " << p1[0] << " " << p1[2];
    out << " L " << p2[0] << " " << p2[2];
  }
  out << "\" />" << endl;
}

void SVGWriter::WriteTree(ostream& out) {

  // draw creases in dotted grey as a single path
  out << "  <path class=\"t\" d=\"";

  for (int i = 0; i < model_->e_size; ++i) {
    const edge& e = model_->edges[i];
    if (e.type == 'b')
      continue;

    // Center of faces
    Vector3d c[2];

    for (int i = 0; i < 2; ++i) {
      const triangle& face = model_->tris[e.fid[i]];
      for (int j = 0; j < 3; ++j)
        c[i] += Vector3d(model_->vertices[face.v[j]].p.get());
      c[i] *= (1.0 / 3);
    }

    Vector3d p1 = GetSVGCoord(c[0]);
    Vector3d p2 = GetSVGCoord(c[1]);

    out << " M " << p1[0] << " " << p1[2];
    out << " L " << p2[0] << " " << p2[2];

  }

  out << "\" />" << endl;
}

void SVGWriter::WriteCreasesInColor(ostream& out) {

  for (int i = 0; i < model_->e_size; ++i) {
    const edge& e = model_->edges[i];
    if (e.type == 'b' || fabs(e.folding_angle) <= 1e-3)
      continue;

    string classname =
        e.folding_angle > 0 ? "m" : (e.folding_angle < 0 ? "v" : "n");

    out << "  <path class=\"" << classname << "\" d=\"";

    Vector3d p1 = GetSVGCoord(e.vid[0]);
    Vector3d p2 = GetSVGCoord(e.vid[1]);

    out << " M " << p1[0] << " " << p1[2];
    out << " L " << p2[0] << " " << p2[2];

    out << "\" />" << endl;
  }

}

void SVGWriter::WriteFlatCreases(ostream& out) {
  for (int i = 0; i < model_->e_size; ++i) {
    const edge& e = model_->edges[i];
    if (e.type == 'b' || fabs(e.folding_angle) > 1e-3)
      continue;

    out << "  <path class=\"n\" d=\"";

    Vector3d p1 = GetSVGCoord(e.vid[0]);
    Vector3d p2 = GetSVGCoord(e.vid[1]);

    out << " M " << p1[0] << " " << p1[2];
    out << " L " << p2[0] << " " << p2[2];

    out << "\" />" << endl;
  }
}

void SVGWriter::AnalyzeNet() {
  this->border_cut_eids_.clear();
  this->zip_line_vids_.clear();
  this->zip_line_eids_.clear();

  //analyze net
  //setup net analyzer
  NetAnalyzer netanalyzer(model_, config_);
  {
    vector<Point2d> boundary_points; //collect boundary vertex positions
    for (int i : boundary_vertices_) {
      Vector3d v = this->GetSVGCoord(i);
      boundary_points.push_back(Point2d(v[0], v[2]));
    }
    //lastone is repeated, so discard
    boundary_points.pop_back();
    netanalyzer.analyze(this->boundary_vertices_, boundary_points);
  }

  //find all border cut edges
  auto border_eids = netanalyzer.getBorderCutEdges();
  this->border_cut_eids_ = border_eids;
  for (int i = 0; i < model_->e_size; ++i) {
    if (model_->edges[i].parent_id == UINT_MAX)
      continue;
    if (find(border_eids.begin(), border_eids.end(), model_->edges[i].parent_id)
        != border_eids.end()) {
      this->border_cut_eids_.push_back(i);
    }
  } //end for i

  //remember the easy start vertices
  this->zip_line_vids_ = netanalyzer.getZipVertices();
  this->zip_line_eids_ = netanalyzer.getZipEdges();

  //copy polygon
  //this->boundary_polygon_.destroy();
  //this->boundary_polygon_.copy(netanalyzer.getNetPolygon());
}

void SVGWriter::WriteZipHints(ostream& out) {
  //
  const float boxw = width_pixel_ * 0.0025;

  out << "  <path class=\"z\" fill=\"none\" d=\"";

  for (list<uint>& zipline : this->zip_line_vids_) {
    for (uint vid : zipline) {
      Vector3d p = GetSVGCoord(vid);
      //draw a small star around each vertex
      if (vid == zipline.front())
        out << " M " << p[0] << " " << p[2];
      else
        out << " L " << p[0] << " " << p[2];
    }
  }
  out << "\" />" << endl;

  //
  //draw a square around the zip start
  out << "  <path class=\"z\" d=\"";
  for (list<uint>& zipline : this->zip_line_vids_) {
    int mid = zipline.size() / 2;
    auto ptr = zipline.begin();
    for (int i = 0; i < mid; i++, ptr++)
      ;
    Vector3d p = GetSVGCoord(*ptr);
    out << " M " << p[0] - boxw << " " << p[2] - boxw;
    out << " L " << p[0] + boxw << " " << p[2] - boxw;
    out << " L " << p[0] + boxw << " " << p[2] + boxw;
    out << " L " << p[0] - boxw << " " << p[2] + boxw;
    out << " L " << p[0] - boxw << " " << p[2] - boxw;
  }
  out << "\" />" << endl;
}

//extra hints for connecting cut edges
void SVGWriter::WriteCutEdgeHints(ostream& out) {
  //hint segments
  list<pair<Vector3d, Vector3d> > hints;

  //
  const float hint_length = width_pixel_ * 0.005;

  for (int i = 0; i < model_->e_size; ++i) {
    const edge& e = model_->edges[i];
    if (e.type != 'b')
      continue; //not a cut edge
    if (e.parent_id == UINT_MAX)
      continue; //no parent ID defined
    if (find(zip_line_eids_.begin(), zip_line_eids_.end(), e.parent_id)
        != zip_line_eids_.end())
      continue; //this edge is a zip line edge
    const edge& pe = model_->edges[e.parent_id]; //parent edge
    assert(pe.type == 'b');           //make sure that we have the right parent
    assert(pe.parent_id == UINT_MAX); //make sure that we have the right parent

    const triangle& t1 = model_->tris[e.fid.front()];
    Vector3d p1 = GetSVGCoord(e.vid[0]);
    Vector3d p2 = GetSVGCoord(e.vid[1]);
    Vector3d p3 = GetSVGCoord(otherv(t1, e));

    const triangle& t2 = model_->tris[pe.fid.front()];
    Vector3d p4 = GetSVGCoord(pe.vid[0]);
    Vector3d p5 = GetSVGCoord(pe.vid[1]);
    Vector3d p6 = GetSVGCoord(otherv(t2, pe));

    //ensure the order is correct...
    auto& e_vid1 = model_->vertices[e.vid[0]];
    auto& e_vid2 = model_->vertices[e.vid[1]];
    auto& pe_vid1 = model_->vertices[pe.vid[0]];
    auto& pe_vid2 = model_->vertices[pe.vid[1]];

    uint evid1_id = getrootid(model_, e.vid[0]); //(e_vid1.parent_id==UINT_MAX)?e.vid[0]:e_vid1.parent_id;
    uint evid2_id = getrootid(model_, e.vid[1]); //(e_vid2.parent_id==UINT_MAX)?e.vid[1]:e_vid2.parent_id;
    uint pevid1_id = getrootid(model_, pe.vid[0]); //(pe_vid1.parent_id==UINT_MAX)?pe.vid[0]:pe_vid1.parent_id;
    uint pevid2_id = getrootid(model_, pe.vid[1]); //(pe_vid2.parent_id==UINT_MAX)?pe.vid[1]:pe_vid2.parent_id;

    if (evid1_id == pevid2_id && evid2_id == pevid1_id) {
      swap(p4, p5);
    } else if (evid1_id == pevid1_id && evid2_id == pevid2_id) {
      //do nothing
    } else {
      cerr << "! Error: Something is not right: (" << evid1_id << ","
          << evid2_id << ") vs. (" << pevid1_id << "," << pevid2_id << ")"
          << endl;
      continue;
    }

    //compute hints
    float s = drand48() * 0.8 + 0.1; //0.1~0.9
    Vector3d m1 = p1 * (1 - s) + p2 * s;
    Vector3d m2 = p4 * (1 - s) + p5 * s;

#if 1
    //cout<<"S="<<s<<endl;
    Vector3d dir1 = (p3 - m1);
    if (dir1.norm() / 2 < hint_length)
      continue; //too short
    dir1 = dir1.normalize() * hint_length;

    Vector3d dir2 = (p6 - m2);
    if (dir2.norm() / 2 < hint_length)
      continue; //too short
    dir2 = dir2.normalize() * hint_length;
#else
    float t=drand48()*0.8+0.1; //0.1~0.9
    if(s>t) swap(s,t);

#endif

    //add hints
    hints.push_back(make_pair(m1, m1 + dir1));
    hints.push_back(make_pair(m2, m2 + dir2));
  }

  sort_line_segments(hints);

  out << "  <path class=\"h\" d=\"";
  for (auto& hint : hints) {
    const auto& p1 = hint.first;
    const auto& p2 = hint.second;
    out << " M " << p1[0] << " " << p1[2];
    out << " L " << p2[0] << " " << p2[2];
  }
  out << "\" />" << endl;
}

void SVGWriter::WriteLabels(ostream & out) {
  for (int i = 0; i < this->model_->e_size; ++i) {
    const edge& e = this->model_->edges[i];
    if (e.type != 'b')
      continue; //not a cut edge

    int id = (e.parent_id == UINT_MAX) ? i : e.parent_id; //find the original id of this edge
    const triangle& t = model_->tris[e.fid.front()];

    Vector3d p0 = GetSVGCoord(e.vid[0]);
    Vector3d p1 = GetSVGCoord(e.vid[1]);
    Vector3d p2 = GetSVGCoord(otherv(t, e));

    auto vec = (p1 - p0).normalize();
    if ((vec % (p2 - p0).normalize())[1] < 0) {
      vec = -vec;
      swap(p0, p1);
    }
    const auto angle = RadToDeg(atan2(vec[2], vec[0])); // angle from x axis in degree

    const auto edge_length = (p1 - p0).norm();
    int label_size = std::to_string(id).size();
    double label_length = label_size * this->font_size_ * 0.4;
    double pp = min(label_length / 2 / edge_length, 0.5);
    const auto center = p0 + (p1 - p0) * (0.5 - pp);

    stringstream ss;
    ss << "<text x=\"" << center[0] << "\" y=\"" << center[2]
        << "\" transform=\"rotate(" << angle << " " << center[0] << " "
        << center[2]
        << ")\" fill=\"darkgreen\" font-weight=\"bold\" font-size=\""
        << this->font_size_ << "\">" << id << "</text>";
    //  sprintf(text_buf,
    //      "<text x=\"%f\" y=\"%f\" transform=\"rotate(%f %f %f)\" fill=\"darkgreen\" font-weight=\"bold\" font-size=\"%f\">%d</text>",
    //      center[0], center[2],   // x, y
    //      angle, center[0], center[2], font_size, id);

    out << ss.str() << endl;

  }  //end for i
}

//assume v1 and v2 are normalized
#if 0 //not used
inline Quaternion rotation(Vector3d v1, Vector3d v2)
{
  v1=v1.normalize();
  v2=v2.normalize();
  Vector3d axis;
  ///const Vector3d axis(0,1,0); //(v1%v2).normalize();
  if((v1%v2)[1]>0) axis=Vector3d(0,1,0);
  else axis=Vector3d(0,-1,0);
  const auto angle=acos(v1*v2);
  return Quaternion::get(angle, axis);
}
#endif

bool SVGWriter::IsValid(const SVGWriter::Tab& tab) const {
  //check intersections with all cut edges
  for (int i = 0; i < this->model_->e_size; ++i) {
    const edge& e = this->model_->edges[i];
    if (e.type != 'b')
      continue; //not a cut edge
    if (i == tab.eid)
      continue; //same edge, avoid checking for intersection

    Vector3d u = GetSVGCoord(e.vid[0]);
    Vector3d v = GetSVGCoord(e.vid[1]);
    if (tab.intersect(u, v)) {
      //cout<<"tab for e="<<tab.eid<<" ("<< model_->edges[tab.eid].parent_id<<") collide with eid="<<eid<<" ("<< model_->edges[eid].parent_id<<")"<<endl;
      return false;
    }

    v = u;
  }

  //check intersections with all existing tabs
  uint tsize = this->tabs_.size();
  for (auto& t : tabs_) {
    if (tab.intersect(t.second))
      return false;
  }

  return true; //no intersections
}

SVGWriter::Tab SVGWriter::BuildTab(Vector3d& a, Vector3d& b, Vector3d& dir,
    double len, uint eid) {
  if (dir.norm() < len) {
    dir = dir.normalize() * len;
  }

  Tab tab;
  Bezier(a, (a + b) / 2 + dir, b, tab.shape_, 10);

  tab.eid = eid;

  return tab;
}

void SVGWriter::BuildTabs() {
  //this function build tabs only for hard-to fold boundary cut edges
  //const float tab_length=width_pixel_ * 0.05;

  for (int i = 0; i < this->model_->e_size; ++i) {
    const edge& e = this->model_->edges[i];
    if (e.type != 'b')
      continue; //not a cut edge
    if (e.parent_id == UINT_MAX)
      continue; //no parent ID defined

    const edge& pe = model_->edges[e.parent_id]; //parent edge
    assert(pe.type == 'b');           //make sure that we have the right parent
    assert(pe.parent_id == UINT_MAX); //make sure that we have the right parent

    //use hard-to fold border cut only
    if (find(border_cut_eids_.begin(), border_cut_eids_.end(), i)
        == border_cut_eids_.end())
      continue;

    //get triangles incident to the cut edges
    const triangle& t1 = model_->tris[e.fid.front()];
    Vector3d p1 = GetSVGCoord(e.vid[0]);
    Vector3d p2 = GetSVGCoord(e.vid[1]);
    Vector3d p3 = GetSVGCoord(otherv(t1, e));

    auto v1 = (p2 - p1).normalize();
    if ((v1 % (p3 - p1).normalize())[1] < 0) {
      swap(p1, p2);
      v1 = -v1;
    }

    const triangle& t2 = model_->tris[pe.fid.front()];
    Vector3d p4 = GetSVGCoord(pe.vid[0]);
    Vector3d p5 = GetSVGCoord(pe.vid[1]);
    Vector3d p6 = GetSVGCoord(otherv(t2, pe));

    //find rotation that will bring pe to e
    auto v2 = (p5 - p4).normalize();
    if ((v2 % (p6 - p4).normalize())[1] < 0) {
      swap(p4, p5);
      v2 = -v2;
    }

    //Quaternion rotq=rotation( p2-p1, p4-p5);

    //tab 1 is p1, m1, (p1+m1)/2+v2
    Vector3d m1 = (p1 + p2) / 2; //make a tab between p1 and m1
    Vector3d m2 = (p2 + m1) / 2;
    Vector3d n1(-v1[2], 0, v1[0]);
    //v1=-v1;
    //Vector3d v1  = rotq.rotate( (p3-m2) )/2; //.normalize();
    //Vector3d v1 = ()

    Vector3d m3 = (p4 + p5) / 2; //make a tab between p4 and m3
    Vector3d m4 = (p5 + m3) / 2;
    Vector3d n2(-v2[2], 0, v2[0]);
    //v2=-v2;
    //Vector3d v2  = (-rotq).rotate( (p6-m4) )/2;//.normalize();

    double tab_length = (p1 - m1).norm();
    Tab tab1 = BuildTab(p1, m1, n1, tab_length, i);
    Tab tab2 = BuildTab(p4, m3, n2, tab_length, e.parent_id);

    if (IsValid(tab1))
      tabs_[tab1.eid] = tab1;
    if (IsValid(tab2))
      tabs_[tab2.eid] = tab2;
  }    //end for i
}

//extra hints for connecting cut edges
void SVGWriter::WriteTabs(ostream& out) {
  out << "  <path class=\"tab\" d=\"";

  for (auto& t : tabs_) {
    auto& tab = t.second;
    const auto& start = tab.shape_.front();
    out << " M " << start[0] << " " << start[2];
    for (auto it = ++tab.shape_.begin(); it != tab.shape_.end(); it++) {
      out << " L " << (*it)[0] << " " << (*it)[2];
    }
    out << " L " << start[0] << " " << start[2];
  }
  out << "\" />" << endl;
}

void SVGWriter::WriteFooter(ostream& out) {
  out << "</svg>" << endl;
}

void SVGWriter::FindBoundaryPolygon(vector<int>* boundary_vertices) {

  boundary_vertices->clear();
  set<int> eids;

  int current_eid = -1;

  // Get the start vertex
  for (int eid = 0; eid < model_->e_size; ++eid) {
    const edge& e = model_->edges[eid];
    if (e.type == 'b') {
      current_eid = eid;
      eids.insert(eid);
      boundary_vertices->push_back(e.vid[0]);
      boundary_vertices->push_back(e.vid[1]);
      break;
    }
  }

  if (current_eid < 0) {
    cerr << "!Error! Cannot find boundary edges in the model" << endl;
    exit(1);
  }

  while (boundary_vertices->front() != boundary_vertices->back()) {
    const vertex& v = model_->vertices[boundary_vertices->back()];
    for (int eid : v.m_e) {
      const edge& e = model_->edges[eid];
      if (e.type != 'b')
        continue;
      if (eids.count(eid))
        continue;

      int vid0 = boundary_vertices->back();
      int vid1 = e.vid[0] == vid0 ? e.vid[1] : e.vid[0];

      // TODO(zxi) fix bug for ghost edges in model::cutEdge
      bool has_edge = false;
      for (int fid : e.fid) {
        const triangle& f = model_->tris[fid];
        set<uint> vs(f.v, f.v + 3);
        if (vs.count(vid0) && vs.count(vid1)) {
          has_edge = true;
          break;
        }
      }

      if (!has_edge) {
        cerr << "!Error! No edges when finding boundary of the polygon!"
            << endl;
        continue;
      }

      boundary_vertices->push_back(vid1);
      eids.insert(eid);
      break;
    }
  }
}

Vector3d SVGWriter::GetSVGCoord(uint vid) const {
  return GetSVGCoord(this->model_->vertices[vid].p);
}

Vector3d SVGWriter::GetSVGCoord(const Vector3d& v) const {
  return GetSVGCoord(Point3d(v.get()));
}

Vector3d SVGWriter::GetSVGCoord(const Point3d& v) const {
  return (v - net_min_) * config_.scale + PADDING;
}

inline double getPt(double n1, double n2, float perc) {
  double diff = n2 - n1;

  return n1 + (diff * perc);
}

void Bezier(const Vector3d & a, const Vector3d & b, const Vector3d & c,
    list<Vector3d>& curve, int size) {
  auto x1 = a[0];
  auto y1 = a[2];
  auto x2 = b[0];
  auto y2 = b[2];
  auto x3 = c[0];
  auto y3 = c[2];

  curve.push_back(a);

  for (int i = 1; i < size; i++) {
    double s = i * 1.0 / size;

    auto xa = getPt(x1, x2, s);
    auto ya = getPt(y1, y2, s);
    auto xb = getPt(x2, x3, s);
    auto yb = getPt(y2, y3, s);

    auto x = getPt(xa, xb, s);
    auto y = getPt(ya, yb, s);

    curve.push_back(Vector3d(x, 0, y));
  }

  curve.push_back(c);
}

} /* namespace util */
} /* namespace unfolding */
} /* namespace masc */
