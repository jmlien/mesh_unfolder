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
#include <unordered_set>

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
  if (this->config_.svg_add_tabs) {
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
      + std::to_string(stroke_width * 2) + ";fill:none;stroke-linejoin:round";
  const string normal_style = "stroke:rgb(128,128,128);stroke-width:"
      + std::to_string(stroke_width);
  const string mountain_style = "stroke:rgb(255,0,0);stroke-width:"
      + std::to_string(stroke_width) + ";stroke-linecap:round;"; // + dashed_line;

  const string valley_style = "stroke:rgb(0,0,255);stroke-width:"
      + std::to_string(stroke_width) + ";stroke-linecap:round;";

  const string hint_style = "stroke:rgb(240,200,0);stroke-width:"
      + std::to_string(stroke_width * 2) + ";stroke-linecap:round;";

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

  const string chamfer_style = "stroke:rgb(0,125,0);stroke-width:"
      + std::to_string(stroke_width / 2) + ";fill:green;fill-opacity:0.5";

  styles_["m"] = mountain_style;
  styles_["v"] = valley_style;
  styles_["b"] = boundary_style;
  styles_["n"] = normal_style;
  styles_["t"] = tree_style;
  styles_["c"] = crease_style;
  styles_["h"] = hint_style;
  styles_["z"] = zip_style;
  styles_["tab"] = tab_style;
  styles_["chamfer"] = chamfer_style;

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
  if (this->config_.svg_add_tabs)
    BuildTabs();

  if (this->config_.svg_valley_chamfer)
    BuildChamfers();

  //don't write anything else if the texture is the only thing we need
  if(type != ExportSVGType::TEXTURE)
  {
    this->WriteBoundaryPolygon(out);

    //chamfers
    if(type & ExportSVGType::CHAMFER && this->config_.svg_valley_chamfer)
    {
      if (this->config_.svg_valley_chamfer)
        WriteChamfers(out);
    }

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
        if (this->config_.svg_add_tabs)
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
  out << "<svg width=\"" << width_pixel_ << "mm\" height=\"" << height_pixel_
      << "mm\"" << " viewBox=\"0 0 "<<width_pixel_<<" "<< height_pixel_<<"\""
      <<" version=\"1.1\" xmlns=\"http://www.w3.org/2000/svg\" "
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
    //if (i != 0 && tabs_.empty() == false) {
    if (tabs_.empty() == false)
    {
      int j = (i + 1) % vsize;
      uint eid = this->geteid(boundary_vertices_[i], boundary_vertices_[j]);
      //find if there is a tab created for this edge
      if (tabs_.find(eid) != tabs_.end())
      {
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
    {
      if(HasTab(i)==false) continue;
    }

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
    if (model_->edges[i].cut_twin_id == UINT_MAX)
      continue;
    if (find(border_eids.begin(), border_eids.end(), model_->edges[i].cut_twin_id)
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
  unordered_set<uint> mids;
  for (list<uint>& zipline : this->zip_line_vids_) {
    int mid = zipline.size() / 2;
    auto ptr = zipline.begin();
    for (int i = 0; i < mid; i++, ptr++);
    mids.insert(*ptr);
  }

  list<Vector3d> mid_pts;
  for(uint mid : mids)
  {
    mid_pts.push_back(GetSVGCoord(mid));
  }
  sort_points(mid_pts);

  out << "  <path class=\"z\" d=\"";
  for(Vector3d & p : mid_pts)
  {
    WriteSquare(out,p,boxw);
    // out << " M " << p[0] - boxw << " " << p[2] - boxw;
    // out << " L " << p[0] + boxw << " " << p[2] - boxw;
    // out << " L " << p[0] + boxw << " " << p[2] + boxw;
    // out << " L " << p[0] - boxw << " " << p[2] + boxw;
    // out << " L " << p[0] - boxw << " " << p[2] - boxw;
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
    if (e.cut_twin_id == UINT_MAX)
      continue; //no parent ID defined
    if (find(zip_line_eids_.begin(), zip_line_eids_.end(), e.cut_twin_id)
        != zip_line_eids_.end())
      continue; //this edge is a zip line edge
    const edge& pe = model_->edges[e.cut_twin_id]; //parent edge
    assert(pe.type == 'b');           //make sure that we have the right parent
    assert(pe.cut_twin_id == UINT_MAX); //make sure that we have the right parent

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

    uint evid1_id  = getrootvid(model_, e.vid[0]);
    uint evid2_id  = getrootvid(model_, e.vid[1]);
    uint pevid1_id = getrootvid(model_, pe.vid[0]);
    uint pevid2_id = getrootvid(model_, pe.vid[1]);

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

void SVGWriter::WriteLabels(ostream & out)
{
  //face id
  for (int i = 0; i < this->model_->t_size; ++i)
  {
    const triangle& t = this->model_->tris[i];
    int id=(t.source_fid==-1)?i:t.source_fid;
    Vector3d p0 = GetSVGCoord(t.v[0]);
    Vector3d p1 = GetSVGCoord(t.v[1]);
    Vector3d p2 = GetSVGCoord(t.v[2]);
    Vector3d center( (p0[0]+p1[0]+p2[0])/3, 0, (p0[2]+p1[2]+p2[2])/3);
    stringstream ss;
    ss << "<text x=\"" << center[0] << "\" y=\"" << center[2]
        << "\" fill=\"gold\" font-weight=\"bold\" font-size=\""
        << this->font_size_ << "\">" << id << "</text>\n";
    out << ss.str() << endl;
  }

  //edge id
  for (int i = 0; i < this->model_->e_size; ++i) {
    const edge& e = this->model_->edges[i];

#if 0 //change to 1 to prevent printing crease lines
    if (e.type != 'b')
     continue; //not a cut edge
#endif

    //int id = (e.cut_twin_id == UINT_MAX) ? i : e.cut_twin_id;
    int id = (e.source_eid == UINT_MAX) ? i : e.source_eid;

    //find the original id of this edge
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

//check if a tab is collision free
bool SVGWriter::IsValid(const SVGWriter::Tab& tab) const
{
  //check intersections with all cut edges
  for (int i = 0; i < this->model_->e_size; ++i)
  {
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

//chamfer related Methods
void SVGWriter::BuildChamfers()
{

    for (int i = 0; i < this->model_->e_size; ++i)
    {
      const edge& e = this->model_->edges[i];
      //if (e.folding_angle >=0 || e.type == 'b') continue; //not a valley

      if (e.type == 'b' || e.folding_angle >=0 ) continue; //a cut edge

      //get triangles incident to the fold edge
      const triangle& t1 = model_->tris[e.fid.front()];
      Vector3d p1 = GetSVGCoord(e.vid[0]);
      Vector3d p2 = GetSVGCoord(e.vid[1]);
      Vector3d p3 = GetSVGCoord(otherv(t1, e));

      auto v1 = (p2 - p1).normalize();
      if ((v1 % (p3 - p1).normalize())[1] < 0) {
        swap(p1, p2);
        v1 = -v1;
      }

      const triangle& t2 = model_->tris[e.fid.back()];
      Vector3d p4 = GetSVGCoord(otherv(t2, e));

      //so the incident polygon is p1, p4, p2, p3
      chamfers_[i]=BuildChamfer(p1,p4,p2,p3, config_.svg_valley_chamfer_width_ratio, i);
    }    //end for i
}

//build a chamfer around edge (a,c) with opposite vertices c and d
//width of the angle is porpotion to the dihedral angle of eid
SVGWriter::Chamfer SVGWriter::BuildChamfer
(Vector3d& a, Vector3d& b, Vector3d& c, Vector3d& d, float len_ratio, uint eid)
{
  const edge& e = this->model_->edges[eid];
  float width= fabs(e.folding_angle*len_ratio);

  Vector3d ca=(c-a).normalize();
  Vector3d ba=(b-a);
  Vector3d bc=(b-c);
  Vector3d n(ca[2],ca[1],-ca[0]);
  double W_ba=width/fabs(ba*n);
  double W_bc=width/fabs(bc*n);
  if(W_ba>1) W_ba=1;
  if(W_bc>1) W_bc=1;
  Vector3d b1=a+ba*W_ba;
  Vector3d b2=c+bc*W_bc;

  Vector3d da=(d-a);
  Vector3d dc=(d-c);
  double W_da=width/fabs(da*n);
  double W_dc=width/fabs(dc*n);
  if(W_da>1) W_da=1;
  if(W_dc>1) W_dc=1;
  Vector3d d1=a+da*W_da;
  Vector3d d2=c+dc*W_dc;

  //build chamfer
  Chamfer chamfer;
  chamfer.eid=eid;
  chamfer.width=width;
  chamfer.hex_[0]=a;
  chamfer.hex_[1]=b1;
  chamfer.hex_[2]=b2;
  chamfer.hex_[3]=c;
  chamfer.hex_[4]=d2;
  chamfer.hex_[5]=d1;

  return chamfer;
}

//create a tab for a given edge, where a and b defines the edges
//dir is perpendicular to ab and len the
SVGWriter::Tab SVGWriter::BuildTab
(Vector3d& a, Vector3d& b, Vector3d& dir,double len, uint eid)
{
  if (dir.norm() < len) {
    dir = dir.normalize() * len;
  }

  Tab tab;
  Bezier(a, (a + b) / 2 + dir, b, tab.shape_, 3);

  tab.eid = eid;

  return tab;
}

void SVGWriter::BuildTabs() {
  //this function build tabs only for hard-to fold boundary cut edges
  //const float tab_length=width_pixel_ * 0.05;

  for (int i = 0; i < this->model_->e_size; ++i)
  {
    const edge& e = this->model_->edges[i];

    if (e.type != 'b')
      continue; //not a cut edge

    if (e.cut_twin_id == UINT_MAX) //no twin ID defined
    {
      model * src_m = this->model_->source;
      if(src_m==NULL) //no source model
        continue;
      else
      {
        //we may still build a tab if we cannot find the twin edge
        //as this edge might be cut by net surgent, in which this->model is
        //a subset of this->model->source
        const edge & src_e = src_m->edges[e.source_eid];
        const triangle & e_incident_tri = this->model_->tris[e.fid.front()];
        uint src_fid=e_incident_tri.source_fid;
        uint src_other_fid = src_e.otherf(src_fid);
        if(src_other_fid==UINT_MAX ) continue; //no other f... so e is not a cut edge
        if(src_other_fid>src_fid ) continue; //avoid building tabs on both edges

        //build the tab around this edge
        //get triangles incident to the cut edges
        const triangle& t1 = model_->tris[e.fid.front()];
        Vector3d p1 = GetSVGCoord(e.vid[0]);
        Vector3d p2 = GetSVGCoord(e.vid[1]);
        Vector3d p3 = GetSVGCoord(otherv(t1, e));
        //the vector along e
        auto v1 = (p2 - p1).normalize();
        if ((v1 % (p3 - p1).normalize())[1] < 0) {
          swap(p1, p2);
          v1 = -v1;
        }
        //find rotation that will bring pe (i.e, v2) to e (i.e., v1)
        //tab 1 is p1, m1, (p1+m1)/2+v2

        Vector3d m1 = (p1 + p2) / 2; //make a tab between p1 and m1
        //Vector3d m2 = (p2 + m1) / 2;

        Vector3d n1(-v1[2], 0, v1[0]);
        double elen=(m1-p1).norm();
        //double tab_length1 = min(elen,fabs((p6 - p4)*n2.normalize()*0.9));
        double tab_length = min(elen,fabs((p3 - p1)*n1.normalize()*0.9));
        Tab tab = BuildTab(p1, p2, n1, tab_length, i);
        tab.inter_net_tab=true; //mark this tab as a special tab

        if (IsValid(tab))
        {
          tabs_[tab.eid] = tab;
        }

        continue;
      }
    }

    const edge& pe = model_->edges[e.cut_twin_id]; //cut twin edge
    assert(pe.type == 'b');           //make sure that we have the right edge
    assert(pe.cut_twin_id == UINT_MAX); //make sure that we have the right edge

    //get triangles incident to the cut edges
    const triangle& t1 = model_->tris[e.fid.front()];
    Vector3d p1 = GetSVGCoord(e.vid[0]);
    Vector3d p2 = GetSVGCoord(e.vid[1]);
    Vector3d p3 = GetSVGCoord(otherv(t1, e));

    //this is the opposite triangle, t2
    const triangle& t2 = model_->tris[pe.fid.front()];
    Vector3d p4 = GetSVGCoord(pe.vid[0]);
    Vector3d p5 = GetSVGCoord(pe.vid[1]);
    Vector3d p6 = GetSVGCoord(otherv(t2, pe));

    //the vector along e
    auto v1 = (p2 - p1).normalize();
    if ((v1 % (p3 - p1).normalize())[1] < 0) {
      swap(p1, p2);
      v1 = -v1;
    }

    //the vector along pe
    auto v2 = (p5 - p4).normalize();
    if ((v2 % (p6 - p4).normalize())[1] < 0) {
      swap(p4, p5);
      v2 = -v2;
    }

    //find rotation that will bring pe (i.e, v2) to e (i.e., v1)
    //tab 1 is p1, m1, (p1+m1)/2+v2
    Vector3d m1 = (p1 + p2) / 2; //make a tab between p1 and m1
    Vector3d m2 = (p2 + m1) / 2;
    Vector3d n1(-v1[2], 0, v1[0]);

    Vector3d m3 = (p4 + p5) / 2; //make a tab between p4 and m3
    Vector3d m4 = (p5 + m3) / 2;
    Vector3d n2(-v2[2], 0, v2[0]);

    //build tabs for hard-to fold border cut
    //if (find(border_cut_eids_.begin(), border_cut_eids_.end(), i) != border_cut_eids_.end())


    double elen=(m1-p1).norm();
    double tab_length1 = min(elen,fabs((p6 - p4)*n2.normalize()*0.9));
    double tab_length2 = min(elen,fabs((p3 - p1)*n1.normalize()*0.9));
    Tab tab1 = BuildTab(p1, p2, n1, tab_length1, i);
    Tab tab2 = BuildTab(p4, p5, n2, tab_length2, e.cut_twin_id); //pe id
    if(tab_length2>tab_length1){ Tab tmp=tab2; tab2=tab1; tab1=tmp; } //swap

    if (IsValid(tab1))
    {
      tabs_[tab1.eid] = tab1;
    }
    else if (IsValid(tab2))
    {
      tabs_[tab2.eid] = tab2;
    }
    else //still failed...
    {
      //check if e and pe share a vertex...the tab will be different if they do
      uint incident_vid=UINT_MAX;
      if(e.vid[0]==pe.vid[0] || e.vid[0]==pe.vid[1]) incident_vid=e.vid[0];
      else if(e.vid[1]==pe.vid[0] || e.vid[1]==pe.vid[1]) incident_vid=e.vid[1];

      if(incident_vid!=UINT_MAX)
      {
        Vector3d px=(pe.vid[0]==incident_vid)?p5:p4; //the triangle p1,p2,px is the empty area
        Tab tab;
        tab.shape_.push_back(p1);
        //tab.shape_.push_back(p1+(px-p1)*0.3+(p2-p1)*0.3);
        tab.shape_.push_back(p2+(px-p2)*0.8);
        tab.shape_.push_back(p2);
        tab.eid=i;
        if (IsValid(tab)) tabs_[tab.eid] = tab;
      }
    }
  }//end for i
}

//extra hints for connecting cut edges
void SVGWriter::WriteTabs(ostream& out)
{
  const float boxw = width_pixel_ * 0.0025;
  out << "  <path class=\"tab\" d=\"";

  for (auto& t : tabs_) {
    auto& tab = t.second;
    const auto& start = tab.shape_.front();
    out << " M " << start[0] << " " << start[2];
    for (auto it = ++tab.shape_.begin(); it != tab.shape_.end(); it++) {
      out << " L " << (*it)[0] << " " << (*it)[2];
    }
  }
  out << "\" />" << endl;

  //compute tab center and draw a cross?
  list<Vector3d> tcs;
  list<Vector3d> tcs_inter;
  for (auto& t : tabs_) {
    auto& tab = t.second;
    Vector3d tc;
    for (auto& pt : tab.shape_) tc=tc+pt;
    tc=tc/tab.shape_.size();
    if(tab.inter_net_tab) tcs_inter.push_back(tc);
    else tcs.push_back(tc);
  }

  sort_points(tcs_inter);
  sort_points(tcs);

  out << "  <path class=\"tab\" d=\"";
  for(auto& tc:tcs)
  {
    //draw a cross around the tab center
    WriteCross(out,tc,boxw);
  }
  out << "\" />" << endl;
  out << "  <path class=\"tab\" d=\"";
  for(auto& tc:tcs_inter)
  {
    //draw a cross around the tab center
    WriteCross(out,tc,boxw);
    //add a square around the cross
    WriteSquare(out,tc,boxw);
  }
  out << "\" />" << endl;
}


//output chamfer info for valley folds
void SVGWriter::WriteChamfers(ostream& out) {

  for (auto& c : chamfers_) {
    auto& chamfer = c.second;
    const auto& start = chamfer.hex_[0];
    out << "  <path class=\"chamfer\" d=\"";
    out << " M " << start[0] << " " << start[2];
    for (short i=1;i<6;i++) {
      out << " L " << chamfer.hex_[i][0] << " " << chamfer.hex_[i][2];
    }
    out << " L " << start[0] << " " << start[2];
    out << "\" />" << endl;
  }

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

//greedy approach to sort points so the total travel distance is minimized
void SVGWriter::sort_points(list< Vector3d >& pts) const
{
  if(pts.empty()) return;
  list< Vector3d > sorted;
  sorted.push_back(pts.front());
  pts.pop_front();

  while(pts.empty()==false)
  {
    //find a closest pt to the last pt in sorted
    auto& last = sorted.back();
    float min_d=FLT_MAX;
    auto best=pts.end();

    for(auto i=pts.begin();i!=pts.end();i++)
    {
       float d=(*i - last).normsqr();
       if(d<min_d)
       {
         min_d=d;
         best=i;
       }
    }//end for

    assert(best!=pts.end());
    sorted.push_back(*best);
    pts.erase(best);

  }//end while

  pts.swap(sorted);
}

//greedy approach to sort line segments so the total travel distance is minimized
void SVGWriter::sort_line_segments(list< pair<Vector3d,Vector3d> >& segs) const
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

//drawing Methods
void SVGWriter::WriteSquare(ostream& out, const Vector3d& tc, float boxw) const
{
  out << " M " << tc[0] - boxw << " " << tc[2] - boxw;
  out << " L " << tc[0] + boxw << " " << tc[2] - boxw;
  out << " L " << tc[0] + boxw << " " << tc[2] + boxw;
  out << " L " << tc[0] - boxw << " " << tc[2] + boxw;
  out << " L " << tc[0] - boxw << " " << tc[2] - boxw;
}

void SVGWriter::WriteCross(ostream& out, const Vector3d& pt, float boxw) const
{
  out << " M " << pt[0] - boxw << " " << pt[2] - boxw;
  out << " L " << pt[0] + boxw << " " << pt[2] + boxw;
  out << " M " << pt[0] - boxw << " " << pt[2] + boxw;
  out << " L " << pt[0] + boxw << " " << pt[2] - boxw;
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

//
//TAB
//

SVGWriter::Tab::Tab()
{
  //setup default values
  this->eid=UINT_MAX;         //eid of the edge that this tab is appended to
  this->inter_net_tab=false; //this tab connects two nets
}

bool SVGWriter::Tab::intersect(const Tab& other) const
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

bool SVGWriter::Tab::intersect(const Vector3d & u, const Vector3d & v) const
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

} /* namespace util */
} /* namespace unfolding */
} /* namespace masc */
