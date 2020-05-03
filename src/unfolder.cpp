/*
 * unfolder.cpp
 *
 *  Created on: Nov 14, 2014
 *      Author: zhonghua
 */

#include "unfolder.h"

#include <cfloat>
#include <climits>
#include <ctime>
#include <queue>
#include <fstream>
#include <map>
#include <numeric>
#include <functional>
#include <algorithm>
#include <string>
using namespace std;

#include "intersection.h"
#include "CD.h"
#include "OverlappingChecker.h"
#include "Splitter.h"
#include "mathtool/Box.h"
#include "mathtool/tri-tri.h"
#include "util/Statistics.h"
#include "util/ConvexHull2D.h"
#include "util/DisjointSets.h"
#include "MeshCutter.h"
using namespace masc;
using namespace masc::util;
using namespace masc::unfolding;

#if WIN32
const float EdgeWeight::DiagnalEdge = -0.01;
const float EdgeWeight::FlatEdge = -0.0001;
const float EdgeWeight::KeepEdge = 0.0;
const float EdgeWeight::CutEdge = 1.0;
#endif

Unfolder::Unfolder(model* m, const Config& config)
{
  this->m_m = m;
  this->m_config = config;
  this->m_base_face_id = UINT_MAX;
  this->m_rotation_angle = 0.0;
  this->m_flat_edges = 0;
  this->m_color = Vector3d(0.8, 0.8, 0.8);
  this->m_is_flattened = false;
  this->m_max_edge_length = FLT_MIN;
  this->m_min_edge_length = FLT_MAX;
  this->m_top_vertex_id = INT_MIN;
  this->m_spliiter = nullptr;
  this->m_cd = new RAPID_CD(this);
  assert(this->m_cd);
  this->m_selected_edges=new bool[this->m_m->e_size];
  assert(this->m_selected_edges);

  this->m_max_path_len = -1;
  this->m_avg_path_len = -1;

  this->m_check_collision_calls = 0;
  this->m_check_overlapping_calls = 0;

  this->m_total_check_collison_time = 0;
  this->m_total_check_overlapping_time = 0;

  this->m_last_overlap_count = 0;

  this->m_cluster_id = 0;

  if (!config.texture_path.empty()) {
    // override texture path
    m->texture_path = config.texture_path;
  }

  if (!m->texture_path.empty()) {
    this->m_texture_renderer.reset(new TextureRenderer2D(m->texture_path));
  }
}

Unfolder::~Unfolder() {
  //delete m_m;//don't delete this, this is from the client
  delete [] m_selected_edges;
  delete m_cd;
}

void Unfolder::measureModel() {
  int count = 0;

  for (int i = 0; i < this->m_m->v_size; ++i)
    if (this->m_m->vertices[i].hyperbolic)
      ++count;

  cout << " - hyperbolic vertices = " << count << " ( "
      << (count * 100.0 / this->m_m->v_size) << " % )" << endl;

  count = 0;
  for (int i = 0; i < this->m_m->e_size; ++i)
    if (this->m_m->edges[i].type == 'r')
      ++count;

  cout << " - concave edges = " << count<<"/"<<this->m_m->e_size << " ( "
      << (count * 100.0 / this->m_m->e_size) << " % )" << endl;

}

void Unfolder::buildUnfolding() {
  cerr << "- unfolding..." << endl;

  auto start = clock();

  // create a splitter
  this->m_spliiter.reset(Splitter::createSplitter(m_config.heuristic));
  this->m_spliiter->measure(this->m_m);

  int max_tries = m_config.max_retries;

  vector<float> best_weights;

  Statistics<float> all_cut_lengths;
  Statistics<int> all_overlaps;
  Statistics<double> all_convex_hull_areas;

  int flat_count = 0;

  for (int i = 0; i < max_tries; i++) {
    auto weights = this->m_spliiter->assignWeights(this->m_m, this->m_config);
    auto count = this->buildFromWeights(weights);

    cerr << "- iter = " << i << ", total overlaps = " << count << "\r"<<flush;

    if (count < all_overlaps.Min()) {
      best_weights = weights;
    }

    all_overlaps.Add(count);

    if (count == 0) {
      cerr << "- no overlapping found in " << (i + 1) << " tries" << endl;
      cerr << "- base face = " << this->m_base_face_id << endl;
      this->m_is_flattened = true;

      auto cut_length = this->getTotalCutLength();

      auto hull = this->buildConvexHull2D();
      double hull_area = masc::util::ConvexHull2DArea(hull);
      cerr << "- total cut length = " << cut_length << " hull area = "
           << hull_area << endl;

      if (cut_length < all_cut_lengths.Min()) {
        best_weights = weights;
      }

      all_cut_lengths.Add(cut_length);
      all_convex_hull_areas.Add(hull_area);

      if (this->m_config.early_stop)
        break;
    }
  }//end for i
  cerr<<endl;//done

  if (!best_weights.empty()) {
    this->buildFromWeights(best_weights);

    float cut_length = this->getTotalCutLength();
    all_cut_lengths.Add(cut_length);
    if (!m_config.quite)
      cerr << "- best cut length = " << all_cut_lengths.Min() << endl;
  }

  // redo everything....

  this->initUnfold();

  this->alignModel();

  this->computeUnfolding();

  this->rebuildModel();

  if (!m_config.quite) {
    cout << "- total time " << (float) (clock() - start) / CLOCKS_PER_SEC
        << " s" << endl;

    cout << " - min/avg/max overlaps = " << all_overlaps.Min() << "/"
        << all_overlaps.Mean() << "/" << all_overlaps.Max() << endl;

    cout << " - min/avg/max cut length = " << all_cut_lengths.Min() << "/"
        << all_cut_lengths.Mean() << "/" << all_cut_lengths.Max() << endl;

    cout << " - min/avg/max hull area = " << all_convex_hull_areas.Min() << "/"
        << all_convex_hull_areas.Mean() << "/" << all_convex_hull_areas.Max()
        << endl;

    cout << " - total fat edges = " << m_flat_edges << endl;
  }
}

int Unfolder::buildFromWeights(const string& path) {

  if (!m_config.quite)
    cout << " - Building Unfolding From Weights File = " << path << endl;

  this->m_weights.clear();

  ifstream fin(path, ifstream::in);

  if (!fin.good()) {
    cerr << " ! Error ! Can't open weights file " << path << endl;
    return INT_MAX;
  }

  float weight;

  while (!fin.eof()) {
    fin >> weight;
    m_weights.push_back(weight);
  }

  // last line is blank
  m_weights.pop_back();

  if (m_weights.size() != this->m_m->e_size) {
    cerr << " ! Error ! Weights size is incorrect! Excepted = " << m_m->e_size
        << " Actual = " << m_weights.size() << endl;
    return INT_MAX;
  }

  // build the model from weights
  return this->buildFromWeights(this->m_weights);
}

/// build the unfolding from binary weights and return # overlaps
/// Note: this could be done faster wihtout calling buildFromWeights(const vector<float>& weights...)
/// as we can quickly build MST using bweights
int Unfolder::buildFromWeights(const BitVector& bweights, bool force_check_overlaps)
{
	vector<float> weights;
	for (int i=0;i<bweights.size();i++) weights.push_back(bweights[i]?0:1);
	return buildFromWeights(weights, force_check_overlaps);
}

/// build the unfolding from weights and return # overlaps
int Unfolder::buildFromWeights(const vector<float>& weights,
    bool force_check_overlaps) {

  this->m_weights = weights;

  GRAPH g;
  vector<int> parents;

  auto start = clock();

  this->initUnfold();
  this->alignModel();

  // assign weights on dual graph
  for (auto i = 0; i < this->m_m->e_size; i++) {
    const auto& edge = m_m->edges[i];
    const auto fid1 = edge.fid[0];
    const auto fid2 = edge.fid[1];

    // skip boundary edges...
    if (edge.type == 'b')
      continue;

    // set the weight
    g[fid1][fid2] = g[fid2][fid1] = weights[i];
  }

  this->buildMST(g);
  this->computeUnfolding();

  int count = INT_MAX;

  if (force_check_overlaps) {
    count = this->checkOverlaps();
  }

  if (!m_config.quite)
    cout << "- total time " << (float) (clock() - start) / CLOCKS_PER_SEC
        << " s" << endl;

  return count;
}

/// Build the unfolding from a cut file. each cut is one line in <vid1, vid2> 1-based format.
int Unfolder::buildFromCuts(const string& path) {
  if (!m_config.quite)
    cerr << " - Building Unfolding From Cut File = " << path << endl;

  ifstream fin(path, ifstream::in);

  if (!fin.good()) {
    cerr << " ! Error ! Can't open cut file " << path << endl;
    return INT_MAX;
  }

  int vid1;
  int vid2;

  vector<pair<int, int>> cuts;

  while (fin >> vid1 >> vid2) {
    cuts.push_back(make_pair(vid1 - 1, vid2 - 1));
  }

  // build the model from cuts
  return this->buildFromCuts(cuts);
}

/// Build the unfolding from a set of cuts, each cut is in <vid1, vid2> 0-based format.
int Unfolder::buildFromCuts(const vector<pair<int, int>>& edges,
    bool checkOverlaps) {

  if (!m_config.quite)
    cerr << " - Building Unfolding From Cuts, # cuts = " << edges.size()
         << endl;

  vector<float> weights(this->m_m->e_size, 0);

  for (const auto& e : edges) {
    const int eid = this->m_m->getEdgeId(e.first, e.second);
    // Cut edges have the highest weights.
    weights[eid] = 1.0;
  }

  return this->buildFromWeights(weights, checkOverlaps);
}

vector<double> Unfolder::getFullFoldingAngles(
    const vector<double>& folding_angles) {
  vector<double> full(this->m_m->e_size);

  auto it = folding_angles.begin();
  for (int i = 0; i < full.size(); ++i) {
    if (this->m_fold_edges.count(i)) {
      full[i] = *it;
      ++it;
    }
  }

  return full;
}

void Unfolder::unfoldTo(const vector<double>& folding_angles) {

  this->m_unfolded = this->m_org;

  vector<Matrix4x4> unfolding_matrax(m_m->t_size);

  //int c = 0;
  for (auto fid : m_ordered_face_list) {

    // no need to unfold base_face;
    if (fid == this->m_base_face_id)
      continue;

    // faces of tabs
    if (fid >= this->m_m->t_size)
      break;

    const auto pfid = m_parents[fid];				// get parent face id
    assert(pfid >= 0 && pfid < this->m_m->t_size);

    const auto& f = this->m_m->tris[fid];				  // current_face
    const auto& pf = this->m_m->tris[pfid];				// parent_face
    const auto eid = this->m_shared_edge[fid][pfid];	// shared edge id
    const auto& e = this->m_m->edges[eid];				// shared edge

    this->m_m->edges[eid].cutted = false;

    // rotation axis order matters, find the correct order
    const auto tid = this->isEdgeCCW(pfid, e.vid[0], e.vid[1]) ? 0 : 1;

    const auto& ev1 = this->getUnfoldedVertex(pfid, e.vid[tid]);// get edge vertex1
    const auto& ev2 = this->getUnfoldedVertex(pfid, e.vid[1 - tid]);// get edge vertex2

    const auto folding_angle = folding_angles[eid];

    assert(!std::isnan(folding_angle));

    // compute unfolding matrix
    unfolding_matrax[fid] = Matrix4x4::getRotationMatrix(ev1, ev2,
        folding_angle) * unfolding_matrax[pfid];

    // unfold current face
    for (auto k = 0; k < 3; k++) {
      this->m_unfolded[fid][k].second = (unfolding_matrax[fid]
          * this->m_unfolded[fid][k].second.from3dto4d()).from4dto3d();
    }
  }
}

void Unfolder::linearUnfoldTo(double percentage, const int count) {

  vector<double> folding_angles(m_m->e_size);

  for (auto fid : m_ordered_face_list) {

    // no need to unfold base_face;
    if (fid == this->m_base_face_id)
      continue;

    // faces of tabs
    if (fid >= this->m_m->t_size)
      break;

    const auto pfid = m_parents[fid];				// get parent face id
    assert(pfid >= 0 && pfid < this->m_m->t_size);

    const auto& f = this->m_m->tris[fid];				// current_face
    const auto& pf = this->m_m->tris[pfid];				// parent_face
    const auto eid = this->m_shared_edge[fid][pfid];	// shared edge id
    const auto& e = this->m_m->edges[eid];				// shared edge

    const auto folding_angle = -e.folding_angle * percentage;

    assert(!std::isnan(folding_angle));

    folding_angles[eid] = folding_angle;
  }

  this->unfoldTo(folding_angles);
}

void Unfolder::orderedUnfoldTo(double percentage) {

  int count = this->m_m->t_size * percentage;

  double ppf = (1.0 / this->m_m->t_size);

  double left_percentage = (percentage - ppf * count) / ppf;

  vector<double> folding_angles(m_m->e_size);

  auto face_list = m_ordered_face_list;

  std::reverse(face_list.begin(), face_list.end());

  int c = 0;

  for (auto fid : face_list) {
    ++c;
    // no need to unfold base_face;
    if (fid == this->m_base_face_id)
      continue;

    // faces of tabs
    if (fid >= this->m_m->t_size)
      break;

    const auto pfid = m_parents[fid];				// get parent face id
    assert(pfid >= 0 && pfid < this->m_m->t_size);

    const auto& f = this->m_m->tris[fid];				// current_face
    const auto& pf = this->m_m->tris[pfid];				// parent_face
    const auto eid = this->m_shared_edge[fid][pfid];	// shared edge id
    const auto& e = this->m_m->edges[eid];				// shared edge

    auto folding_angle = -e.folding_angle;

    // fully unfolded expect the last one
    if (c > count)
      folding_angle *= left_percentage;

    assert(!std::isnan(folding_angle));

    folding_angles[eid] = folding_angle;

    if (c > count)
      break;
  }

  this->unfoldTo(folding_angles);

}

//void Unfolder::foldToCompactState() {
//  auto elist = this->m_ordered_crease_list;
//  std::reverse(elist.begin(), elist.end());
//
//  this->m_evolving_seqs.clear();
//
//  vector<double> folding_angles(this->m_m->e_size, 0);
//
//  for (auto eid : elist) {
//    const auto e = this->m_m->edges[eid];
//    // ignore diagonal edges
//    if (e.diagonal)
//      continue;
//    folding_angles[eid] = -e.folding_angle;
//  }
//
//  this->m_evolving_seqs.push_back(folding_angles);
//
//  for (auto eid : elist) {
//    const auto e = this->m_m->edges[eid];
//    // ignore diagonal edges
//    if (e.diagonal)
//      continue;
//    auto fa = folding_angles[eid];
//    const int steps = 30;
//    for (int i = 0; i < steps; ++i) {
//      folding_angles[eid] = fa - (i + 1.0) / steps * PI;
//      this->m_evolving_seqs.push_back(folding_angles);
//    }
//  }
//}

void Unfolder::unfoldTo(double percentage) {
  if (percentage < 0)
    percentage = 0.0f;
  if (percentage > 1)
    percentage = 1.0f;

  if (m_config.ordered_unfolding) {
    this->orderedUnfoldTo(percentage);
  } else {
    this->linearUnfoldTo(percentage);
  }
}

/////////////////////////////////////////////////////////////////////////////

// get vertices of current folded state
vector<Vector3d> Unfolder::getVertices() {
  vector<Vector3d> vs;

  for (auto f : this->m_unfolded)
    for (auto p : f)
      vs.push_back(p.second);

  return vs;
}

void Unfolder::findMostCompactState() {
  if (!m_config.quite)
    cout << "Finding most compact state..." << endl;

  this->unfoldTo(0.0);
  mathtool::Box3d b_model(this->getVertices());

  if (!m_config.quite)
    cout << "Original mesh: dim = " << b_model.getDim() << " sum = "
        << b_model.getSumOfDim() << endl;

  this->unfoldTo(1.0);
  mathtool::Box3d b_net(this->getVertices());

  if (!m_config.quite)
    cout << "Net: dim = " << b_net.getDim() << " sum = " << b_net.getSumOfDim()
        << endl;

  double min_sum_dim = FLT_MAX;

  // std::min(b_net.getSumOfDim(), b_model.getSumOfDim());

  vector<double> best_cfg(this->m_m->e_size);
  vector<double> curr_cfg(this->m_m->e_size);

  for (int i = 0; i < m_config.run; ++i) {
    //cout << "i = " << i << endl;

    for (auto fid : m_ordered_face_list) {

      // no need to unfold base_face;
      if (fid == this->m_base_face_id)
        continue;

      // faces of tabs
      if (fid >= this->m_m->t_size)
        break;

      const auto pfid = m_parents[fid];       // get parent face id
      assert(pfid >= 0 && pfid < this->m_m->t_size);

      const auto& f = this->m_m->tris[fid];       // current_face
      const auto& pf = this->m_m->tris[pfid];       // parent_face
      const auto eid = this->m_shared_edge[fid][pfid];  // shared edge id
      const auto& e = this->m_m->edges[eid];        // shared edge

      auto folding_angle = 0.0;

      switch (rand() % 6) {
      case 0:
        folding_angle = 0.0;
        break;
      case 1:
        folding_angle = -e.folding_angle;
        break;
      case 2:
        folding_angle = -e.folding_angle * mathtool::drand48();
        break;
      case 3:
        folding_angle = mathtool::drand48() * PI;
        if (rand() % 2 == 1)
          folding_angle *= -1;
        break;
      case 4:
        folding_angle = PI;
        break;
      case 5:
        folding_angle = PI * 0.95;
        if (rand() % 2 == 1)
          folding_angle *= -1;
        break;
      }

      folding_angle = (mathtool::drand48() - 0.5) * 2 * PI;

      if (m_config.less_cuts && e.folding_angle == 0)
        folding_angle = 0.0;

      assert(!std::isnan(folding_angle));

      curr_cfg[eid] = folding_angle;
    }

    this->unfoldTo(curr_cfg);
    mathtool::Box3d b_curr(this->getVertices());

    this->checkCollision();
//      continue;

    if (!m_config.quite)
      cout << "Curr: dim = " << b_curr.getDim() << " sum = "
          << b_curr.getSumOfDim() << endl;

    if (b_curr.getSumOfDim() < min_sum_dim) {
      min_sum_dim = b_curr.getSumOfDim();
      best_cfg = curr_cfg;

      if (!m_config.quite)
        cout << "New best: dim = " << b_curr.getDim() << " sum = "
            << b_curr.getSumOfDim() << endl;
    }
  }

  if (!m_config.quite)
    cout << "best: sum = " << min_sum_dim << endl;

  this->unfoldTo(best_cfg);
}

void Unfolder::shrink()
{
  for (int i = 0; i < this->m_m->t_size; ++i) {
    auto& p0 = this->m_unfolded[i][0].second;
    auto& p1 = this->m_unfolded[i][1].second;
    auto& p2 = this->m_unfolded[i][2].second;

    // TODO use circumcenter instead
    auto center = (p0 + p1 + p2) / 3;

    p0 = (p0 - center) * m_config.shrink_factor + center;
    p1 = (p1 - center) * m_config.shrink_factor + center;
    p2 = (p2 - center) * m_config.shrink_factor + center;
  }
}

//return the polygon representing the boundary of the net
masc::polygon::c_ply Unfolder::findBoundaryPolyon()
{
    //SVGWriter writer(m, config);
    if (m_net.v_size == 0) {
      this->rebuildModel();
    }

    if(m_svg_writer.get()==NULL)
      m_svg_writer.reset(new SVGWriter(&m_net, m_config));

    vector<int> boundary;
    m_svg_writer->FindBoundaryPolygon(&boundary);

    //create netshadow
    masc::polygon::c_ply ply(masc::polygon::c_ply::POUT);
    ply.beginPoly();
    for(auto& id : boundary)
    {
      Vector3d pos=m_svg_writer->GetSVGCoord(id);
      ply.addVertex(pos[0],pos[2]);
    }
    ply.endPoly();

    return ply;
}

void Unfolder::dumpObj(const string& path) const {
  cout << "- dumping obj to " << path << endl;

  ofstream out(path);

  vector<Vector3d> vs(m_net.v_size);

  for (int fid = 0; fid < this->m_net.t_size; ++fid) {
    const triangle& f = m_net.tris[fid];
    for (int i = 0; i < 3; ++i) {
      const Vector3d& p = this->m_unfolded[fid][i].second;
      vs[f.v[i]] = p;
    }
  }

  for (int vid = 0; vid < this->m_net.v_size; ++vid) {
    out << "v " << vs[vid][0] << " " << vs[vid][1] << " " << vs[vid][2] << endl;
  }

  for (auto fid = 0; fid < this->m_net.t_size; ++fid) {
    const triangle& f = m_net.tris[fid];
    out << "f " << (f.v[0] + 1) << " " << (f.v[1] + 1) << " " << (f.v[2] + 1)
        << endl;
  }

  out.close();
}

void Unfolder::dumpCreasePattern(const string& path) {

  cout << "- dumping crease pattern to " << path << endl;

  ofstream out(path);

  // 1. fully unfold current model
  this->unfoldTo(1.0);

  // unfolded vertex coordinates
  vector<Vector3d> vs(this->m_m->v_size);

  // fill the coordinates
  for (int fid = 0; fid < this->m_unfolded.size(); ++fid) {
    const triangle& t = this->m_m->tris[fid];
    for (int i = 0; i < 3; ++i)
      vs[t.v[i]] = this->m_unfolded[fid][i].second;
  }

  // dump vertices
  for (auto& v : vs) {
    // swap y and z such that the cp is on xy-plane
    out << "v " << v[0] << " " << v[2] << " " << v[1] << endl;
  }

  // dump faces
  for (int i = 0; i < this->m_m->t_size; ++i) {
    const triangle& t = this->m_m->tris[i];
    out << "f " << (t.v[0] + 1) << " " << (t.v[1] + 1) << " " << (t.v[2] + 1)
        << endl;
  }

  map<long long, vector<int>> symmetry;

  int cid = 0;

  // dump crease lines
  for (int i = 0; i < this->m_m->e_size; ++i) {
    const edge& e = this->m_m->edges[i];
    if (e.type == 'b')
      continue;
    out << "c " << (e.vid[0] + 1) << " " << (e.vid[1] + 1) << " ";

    // assume flat
    if (fabs(e.folding_angle) < 0.005) {
      out << "3 0" << endl;
    } else {
      // non-flat crease lines
      if (e.type == 'c') {
        out << 1 << " " << abs(e.folding_angle) << endl;

        // ignore smaller difference
        long long fa = abs(e.folding_angle) * 1e6;

        symmetry[fa].push_back(cid);
      } else if (e.type == 'r') {
        out << 2 << " " << -abs(e.folding_angle) << endl;

        long long fa = -abs(e.folding_angle) * 1e6;

        symmetry[fa].push_back(cid);
      }
    }

    ++cid;
  }

  // dump boundaries
  for (int i = 0; i < this->m_m->e_size; ++i) {
    const edge& e = this->m_m->edges[i];
    if (e.type != 'b')
      continue;
    out << "c " << (e.vid[0] + 1) << " " << (e.vid[1] + 1) << " 0 0" << endl;
  }

  // dump symmetry
  for (const auto& kv : symmetry) {
    out << "s";
    for (const auto cid : kv.second) {
      out << " " << (cid + 1);
    }
    out << endl;
  }

  out.close();
}

void Unfolder::dumpWrl(const string& path) const {
  cout << "- dumping wrl to " << path << endl;

  ofstream out(path);

  // dump background
  out << "#VRML V2.0 utf8" << std::endl;
  out << "" << std::endl;
  out << "Group {" << std::endl;
  out << "  children [" << std::endl;
  out << "    Background {" << std::endl;
  out << "     skyAngle [ ]" << std::endl;
  out << "     skyColor [ 1.0 1.0 1.0 ]" << std::endl;
  out << "    }" << std::endl;
  out << "  ]" << std::endl;
  out << "}" << std::endl;
  out << std::endl;

  // <cluster_id, max_vids>
  map<int, int> max_vids;
  // <cluster_id, <old_vid, new_vid>>
  map<int, map<int, int>> vids;
  // <cluster, face[i]>: face[i] = {new_v1, new_v2, new_v3}
  map<int, vector<vector<int>>> ffmap;

  // build the mapping
  for (int i = 0; i < this->m_m->t_size; ++i) {
    const auto& t = this->m_m->tris[i];
    int cid = t.cluster_id;
    vector<int> face_vids;
    for (int k = 0; k < 3; k++) {
      // get the new vid for the cluster
      int vid =
          vids[cid].count(t.v[k]) ? vids[cid].at(t.v[k]) : (max_vids[cid]++);
      vids[cid][t.v[k]] = vid;

      face_vids.push_back(vid);
    }

    ffmap[cid].push_back(face_vids);
  }

  for (const auto& kv : vids) {
    int cid = kv.first;
    map<int, int> vmap; // new_id -> old_id;
    for (const auto& vv : kv.second) {
      vmap[vv.second] = vv.first;
    }
    const auto& fs = ffmap[cid];

    Vector3d diffuse(mathtool::drand48() * 0.8 + 0.2,
        mathtool::drand48() * 0.8 + 0.2, mathtool::drand48() * 0.8 + 0.2);
    Vector3d specular(0.1, 0.1, 0.1);
    Vector3d emissive(0.0, 0.0, 0.0);
    double ambient = 0.0;
    double shiniess = 0.05;
    double transparency = 0.0;

    out << std::endl;

    this->writeWrlHeader((int) vmap.size(), (int) fs.size(), out, diffuse,
        specular, emissive, ambient, shiniess, transparency);

    // vertices
    out << "        coord DEF co Coordinate {" << std::endl;
    out << "          point [" << std::endl;
    for (const auto& vv : vmap) {
      const auto& v = this->m_m->vertices[vv.second];
      out << "                  " << v.p[0] << " " << v.p[1] << " " << v.p[2]
          << "," << std::endl;
    }

    out << "         ]" << std::endl;
    out << "         }" << std::endl;

    // faces
    out << "        coordIndex [ " << std::endl;
    for (const auto& f : fs)
      out << "            " << f[0] << ", " << f[1] << ", " << f[2] << ", -1,"
          << std::endl;

    out << "       ]" << std::endl;
    out << "      }" << std::endl;
    out << "    }" << std::endl;
    out << "  ]" << std::endl;
    out << "}" << std::endl;

  }

  out.close();
}

void Unfolder::writeWrlHeader(int nV, int nT, ofstream& out,
    const Vector3d& diffuse, const Vector3d& p, const Vector3d& e,
    double ambient, double shiniess, double transparency) const {
  out << "# Vertices: " << nV << std::endl;
  out << "# Triangles: " << nT << std::endl;
  out << "" << std::endl;
  out << "Group {" << std::endl;
  out << "  children [" << std::endl;
  out << "    Shape {" << std::endl;
  out << "      appearance Appearance {" << std::endl;
  out << "        material Material {" << std::endl;
  out << "          diffuseColor " << diffuse[0] << " " << diffuse[1] << " "
      << diffuse[2] << std::endl;
  out << "          ambientIntensity " << ambient << std::endl;
  out << "          specularColor " << p[0] << " " << p[1] << " " << p[2]
      << std::endl;
  out << "          emissiveColor " << e[0] << " " << e[1] << " " << e[2]
      << std::endl;
  out << "          shininess " << shiniess << std::endl;
  out << "          transparency " << transparency << std::endl;
  out << "        }" << std::endl;
  out << "      }" << std::endl;
  out << "      geometry IndexedFaceSet {" << std::endl;
  out << "        ccw TRUE" << std::endl;
  out << "        solid TRUE" << std::endl;
  out << "        convex FALSE" << std::endl;
}

void Unfolder::dumpWeights(const string& path) const {
  cout << "- dumping weights to " << path << endl;

  ofstream out(path);

  for (const auto w : this->m_weights)
    out << w << endl;

  out.close();
}

void Unfolder::dumpSVG(const std::string& path, ExportSVGType svg_type) {

  if (m_net.v_size == 0)
  {
    this->rebuildModel();
  }

  m_svg_writer.reset(new SVGWriter(&m_net, m_config));

  m_svg_writer->Save(path, svg_type);
}

void Unfolder::dumpOri(const string& path) const {
  cout << "- dumping ori to " << path << endl;

  ofstream out(path);

  // number of vertices
  out << this->m_net.v_size << endl;

  const double tiny = 1e-11;

  // output coordinates of vertices
  for (const auto& v : m_net.vertices) {
    double v0 = fabs(v.p[0]) < tiny ? 0 : v.p[0];
    double v1 = fabs(v.p[1]) < tiny ? 0 : v.p[1];
    double v2 = fabs(v.p[2]) < tiny ? 0 : v.p[2];

    out << v0 << " " << v1 << " " << v2 << endl;
  }

  unordered_map<uint, int> f2c;

  auto cid = -1;
  // output creases
  out << this->m_net.t_size - 1 << endl;
  for (const edge& e : m_net.edges) {

    if (e.type == 'b')
      continue;

    int fid = e.fid[0];
    int pid = e.fid[1];

    if (m_parents[pid] == fid)
      swap(fid, pid);

    int vid0 = e.vid[0];
    int vid1 = e.vid[1];

    ++cid;

    if (fid == m_base_face_id) {
      pid = -1;
      f2c[fid] = -1;
    } else {
      f2c[fid] = cid;

      if (!m_net.isEdgeCCW(pid, vid0, vid1))
        swap(vid0, vid1);
    }

    out << vid0 << " " << vid1 << " " << e.folding_angle << " " << fid << " "
        << pid << endl;
  }

  auto fid = -1;

  // output faces
  out << this->m_net.t_size << endl;
  for (const triangle& f : m_net.tris) {
    ++fid;
    // coordinates

    out << "3 " << f.v[0] << " " << f.v[1] << " " << f.v[2] << " ";
    out << m_parents[fid] << " " << f2c[fid] << endl;
  }

  // base_face_id
  out << m_base_face_id << endl;

  // ordered face list
  for (const auto fid : m_ordered_face_list) {
    out << fid << endl;
  }

  // rotation axis, rotation angle (in rad)
  out << this->m_rotation_axis << -this->m_rotation_angle << endl;

  // output translation
//  out << "0 0 0" << endl;
  out << -this->m_transilation << endl;

  out.close();
}
//
// /// dump unfolded results to svg file to the given path
// void Unfolder::dumpSVGOld(const string& path, const int type) {
//   ofstream out(path);
//
//   auto min_x = DBL_MAX;
//   auto max_x = DBL_MIN;
//   auto min_y = DBL_MAX;
//   auto max_y = DBL_MIN;
//
//   for (const auto& v : this->m_vs) {
//     min_x = min(min_x, v[0]);
//     max_x = max(max_x, v[0]);
//     min_y = min(min_y, v[2]);
//     max_y = max(max_y, v[2]);
//   }
//
//   const auto scale_factor = m_config.scale;
//
//   const auto width = (max_x - min_x);
//   const auto height = (max_y - min_y);
//   const auto aspect = width / height;
//
//   const int width_pixel = width;
//   const int height_pixel = width_pixel / aspect;
//
//   const auto dashed_length = width_pixel * 0.01;
//
//   const auto stroke_width = 1.0 * max(width_pixel, height_pixel) / 800;
//
//   const auto font_size = min(width_pixel, height_pixel) * 0.02
//       * m_config.label_font_scale;
//
//   const string dashed_line = "stroke-dasharray:" + std::to_string(dashed_length)
//       + ", " + std::to_string(dashed_length) + ";";
//
//   const string boundary_style = "stroke:rgb(0,0,0);stroke-width:"
//       + std::to_string(stroke_width * 2) + ";fill:white";
//   const string normal_style = "stroke:rgb(128,128,128);stroke-width:"
//       + std::to_string(stroke_width);
//   const string mountain_style = "stroke:rgb(255,0,0);stroke-width:"
//       + std::to_string(stroke_width) + ";"; // + dashed_line;
//
//   const string valley_style = "stroke:rgb(0,0,255);stroke-width:"
//       + std::to_string(stroke_width) + ";";
//
//   const string tree_style = "stroke:rgb(0,128,0);stroke-width:"
//       + std::to_string(stroke_width) + ";";
//
//   const string crease_style = "stroke:rgb(128,128,128);stroke-width:"
//       + std::to_string(stroke_width) + ";" + dashed_line;
//
//   const string cut_style = "stroke:rgb(0,0,0);stroke-width:"
//       + std::to_string(stroke_width) + ";" + dashed_line;
//   //
// //      stroke-dasharray:"
// //      + std::to_string(dashed_length) + ", " + std::to_string(dashed_length)
// //      + ", " + std::to_string(dashed_length / 5) + ", "
// //      + std::to_string(dashed_length) + ";";
//
//   out << "<?xml version=\"1.0\"?>" << endl;
//   out << "<svg width=\"" << width_pixel << "\" height=\"" << height_pixel
//       << "\"" << " version=\"1.1\" xmlns=\"http://www.w3.org/2000/svg\" "
//       << "xmlns:xlink=\"http://www.w3.org/1999/xlink\">" << endl;
//
//   out << "<defs>" << endl;
//   out << "<style type=\"text/css\"><![CDATA[" << endl;
//
//   out << ".m {" << mountain_style << "}" << endl;
//   out << ".v {" << valley_style << "}" << endl;
//   out << ".b {" << boundary_style << "}" << endl;
//   out << ".n {" << normal_style << "}" << endl;
//   out << ".t {" << tree_style << "}" << endl;
//   out << ".c {" << crease_style << "}" << endl;
//
//   out << "]]></style>" << endl;
//   out << "</defs>" << endl;
//
//   // dump cuts, in two colors
//   if (type & 1) {
//
//     // draw textures first
//     if (this->m_texture_renderer != nullptr && type == 1) {
//
//       this->m_texture_renderer->SetCanvasSize(width_pixel, height_pixel);
//
//       // For each face
//       for (int fid = 0; fid < m_m->t_size; ++fid) {
//         const auto& face = this->m_unfolded[fid];
//
//         const Vector3d& p0_3d = (face[0].second - m_net_min_v) * m_config.scale;
//         const Vector3d& p1_3d = (face[1].second - m_net_min_v) * m_config.scale;
//         const Vector3d& p2_3d = (face[2].second - m_net_min_v) * m_config.scale;
//
//         // Coordinates on SVG canvas
//         const Vector2d p0(p0_3d[0], p0_3d[2]);
//         const Vector2d p1(p1_3d[0], p1_3d[2]);
//         const Vector2d p2(p2_3d[0], p2_3d[2]);
//
//         const triangle& t = this->m_m->tris[fid];
//
//         // UV coordinates
//         const Vector2d& uv0 = this->m_m->texture_pts[t.vt[0]];
//         const Vector2d& uv1 = this->m_m->texture_pts[t.vt[1]];
//         const Vector2d& uv2 = this->m_m->texture_pts[t.vt[2]];
//
//         cerr << "p = " << p0 << "," << p1 << " " << p2 << std::endl;
//         cerr << "uv = " << uv0 << "," << uv1 << " " << uv2 << std::endl;
//
//         this->m_texture_renderer->Render(p0, p1, p2, uv0, uv1, uv2);
//
//       }
//
//       std::string img_url;
//       this->m_texture_renderer->ExportToPngBase64(&img_url);
//       out << "<image xlink:href=\"" << img_url << "\" width=\"" << width_pixel
//           << "px\" height=\"" << height_pixel << "px\" />\n";
//     }
//
//     // draw boundary in black color
//     out << "<path fill=\"none\" style=\"stroke:black;stroke-width:"
//         << std::to_string(stroke_width * 2) << ";\"";
//     out << " d=\"";
//
//     set<pair<int, int>> draw_edges;
//     set<set<int>> processed_edges;
//
//     for (auto& face : this->m_fs) {
//       for (auto j = 1; j <= 3; ++j) {
//         const auto vid0 = face[j - 1];
//         const auto vid1 = face[j % 3];
//         const auto& p0 = this->m_vs[vid0];
//         const auto& p1 = this->m_vs[vid1];
//
//         bool boundary = true;
//
//         auto edge = make_pair(vid0, vid1);
//         auto edge2 = make_pair(vid1, vid0);
//
//         if (!this->m_crease_lines.count(edge)
//             && !this->m_crease_lines.count(edge2)) {
//
//           out << " M " << p0[0] << " " << p0[2];
//           out << " L " << p1[0] << " " << p1[2];
//         }
//       }
//
//       // draw extra cuts for tabs
//       if (type & 8) {
//
//         const int vids[3] = { face[0], face[1], face[2] };
//
//         const Vector3d p[3] = { m_vs[vids[0]], m_vs[vids[1]], m_vs[vids[2]] };
//         const Vector3d e[3] = { p[1] - p[0], p[2] - p[1], p[0] - p[2] };
//
//         const auto shortest_e = std::min(std::min(e[0].norm(), e[1].norm()),
//             e[2].norm());
//
//         const auto x = shortest_e * m_config.extra_cuts_x;
//         const auto y = shortest_e * m_config.extra_cuts_y;
//
//         for (int i = 0; i < 3; ++i) {
//           int vidc = vids[i];
//           int vidn = vids[(i + 1) % 3];
//           set<int> edge = { vidc, vidn };
//
//           if (processed_edges.count(edge))
//             continue;
//           processed_edges.insert(edge);
//
//           const auto ep = e[(i + 2) % 3];
//           const auto ec = e[i];
//           const auto en = e[(i + 1) % 3];
//
//           const auto tp = acos(en.normalize() * ep.normalize());
//           const auto tc = acos(ec.normalize() * ep.normalize());
//           const auto tn = acos(ec.normalize() * en.normalize());
//
//           const auto h = y * sin(tn);
//           const auto b = y * sin(tp) / sin(tn);
//           const auto a = ec.norm() - b - x;
//           const auto z = y * sin(tn) / sin(tc);
//
//           const auto pbx = p[(i + 1) % 3] - (ec.normalize()) * (b + x);
//           const auto pbp = pbx + en.normalize() * y;
//           const auto pz = p[i] - ep.normalize() * z;
//
//           out << " M " << pbx[0] << " " << pbx[2];
//           out << " L " << pbp[0] << " " << pbp[2];
//           out << " L " << pz[0] << " " << pz[2];
//           out << " L " << p[i][0] << " " << p[i][2];
//           out << " L " << pbx[0] << " " << pbx[2];
//
//           auto e01 = make_pair(face[i], face[(i + 1) % 3]);
//           auto e10 = make_pair(face[(i + 1) % 3], face[i]);
//
//           // check whether is boundary edge
//           if (!this->m_crease_lines.count(e01)
//               && !this->m_crease_lines.count(e10))
//             continue;
//
//           // draw mirrow part
//           // height dir: perpendicular to current edge, ignore y...
//           const auto h_dir = Vector3d(-ec[2], 0, ec[0]);
//           const auto pbpp = pbp + h_dir.normalize() * 2 * h;
//           const auto pbzp = pz + h_dir.normalize() * 2 * h;
//
//           out << " L " << pbpp[0] << " " << pbpp[2];
//           out << " L " << pbzp[0] << " " << pbzp[2];
//           out << " L " << p[i][0] << " " << p[i][2];
//
//         }
//       } // end of if (type & 8) extra cuts
//     }
//
//     out << "\" />" << endl;
//
//     {
//       // draw creases in dotted grey as a single path
//       out << "<path fill=\"none\" class=\"c\" d=\"";
//
//       for (auto& c : this->m_cs) {
//         const auto& p1 = this->m_vs[c.vid1];
//         const auto& p2 = this->m_vs[c.vid2];
//
//         // skip flat edges
//         if (c.folding_angle == 0)
//           continue;
//
//         out << " M " << p1[0] << " " << p1[2];
//         out << " L " << p2[0] << " " << p2[2];
//       }
//
//       out << "\" />" << endl;
//     } // end of creases
//
//     // draw spanning tree
//     if (type & 4) {
//       for (auto& c : this->m_cs) {
//         int fid = c.fid;
//         int pid = c.pid;
//         if (pid == -1)
//           continue;
//         const auto& f1 = this->m_fs[fid];
//         const auto& f2 = this->m_fs[pid];
//         Vector3d c1;
//         Vector3d c2;
//         for (int i = 0; i < 3; ++i) {
//           c1 += this->m_vs[f1[i]];
//           c2 += this->m_vs[f2[i]];
//         }
//         c1 = c1 * (1.0 / 3);
//         c2 = c2 * (1.0 / 3);
//
//         char buf[1024];
//
//         sprintf(buf,
//             "<line x1=\"%f\" y1=\"%f\" x2=\"%f\" y2=\"%f\" class=\"t\" />",
//             c1[0], c1[2], c2[0], c2[2]);
//         out << buf << endl;
//       }
//     } // end of spanning tree
//
//   } // end of type 1
//
//   // display mode
//   // 1. draw crease and boundary in their own colors
//   // 2. draw cuts with different thickness
//   if (type == 2) {
//     set<pair<int, int>> draw_edges;
//
//     const double min_line_width = 0.2;
//     const double max_line_width = 3.0;
//
//     char buf[1024];
//
//     // draw crease lines in color
//     for (auto& face : this->m_fs) {
//       for (auto j = 1; j <= 3; ++j) {
//         int vid1 = face[j - 1];
//         int vid2 = face[j % 3];
//         const auto p0 = this->m_vs[vid1];
//         const auto p1 = this->m_vs[vid2];
//
// //        const auto edge = std::make_pair(min(vid1, vid2), max(vid1,vid2));
//         const auto edge1 = std::make_pair(vid1, vid2);
//         const auto edge2 = std::make_pair(vid2, vid1);
//
//         int org_eid = this->m_vemap[edge1];
//         const auto& org_e = this->m_m->edges[org_eid];
//
//         string class_name = "b";
//         bool boundary = true;
//
//         // already drawn
//         if (draw_edges.count(edge1) || draw_edges.count(edge2))
//           continue;
//
//         // const Crease& crease = this->m_cs[m_crease_lines.at(edge1)];
//         if (org_e.folding_angle > 1e-3)
//           class_name = "m";
//         else if (org_e.folding_angle < -1e-3)
//           class_name = "v";
//         else
//           class_name = "n";
//
//         const double line_width = stroke_width
//             * (min_line_width
//                 + fabs(org_e.folding_angle / PI)
//                     * (max_line_width - min_line_width));
//
//         const string style = "stroke-width:" + std::to_string(line_width) + ";";
//
//         boundary = false;
//
//         draw_edges.insert(edge1);
//         draw_edges.insert(edge2);
//
//         sprintf(buf,
//             "<line x1=\"%f\" y1=\"%f\" x2=\"%f\" y2=\"%f\" class=\"%s\" style=\"%s\" />",
//             p0[0], p0[2], p1[0], p1[2], class_name.c_str(), style.c_str());
//         out << buf << endl;
//       }
//     }
//
//     // draw labels
//     if (this->m_config.svg_dump_labels) {
//       char text_buf[1024];
//
//       for (auto i = 0; i < this->m_m->e_size; ++i) {
//         const auto& new_edges = this->m_eemap.at(i);
//
//         // boundary edges
//         if (new_edges.size() == 1u)
//           continue;
//
//         for (const auto& new_edge : new_edges) {
//           if (this->m_crease_lines.count(new_edge))
//             continue;
//
//           const auto p0 = m_vs[new_edge.first];
//           const auto p1 = m_vs[new_edge.second];
//
//           const auto vec = (p1 - p0).normalize();
//           const auto angle = RadToDeg(atan2(vec[2], vec[0])); // angle from x axis in degree
//
//           const auto edge_length = (p1 - p0).norm();
//           int label_size = std::to_string(i).size();
//           double label_length = label_size * font_size * 0.4;
//           double pp = min(label_length / 2 / edge_length, 0.5);
//           const auto center = p0 + (p1 - p0) * (0.5 - pp);
//
//           sprintf(text_buf,
//               "<text x=\"%f\" y=\"%f\" transform=\"rotate(%f %f %f)\" fill=\"darkgreen\" font-weight=\"bold\" font-size=\"%f\">%d</text>",
//               center[0], center[2],   // x, y
//               angle, center[0], center[2], font_size, i);
//
//           out << text_buf << endl;
//
//         } // end for
//
//       } // end for
//
//     } // end if
//   }
//
//   out << "</svg>" << endl;
//
//   out.close();
//
//   cout << "- dumped svg to " << path << endl;
// }

// =============================================================
// private
// =============================================================

void Unfolder::buildDualGraph(GRAPH& g) {
  auto start = clock();

  if (!m_config.quite) {
    cout << "- building dual graph...";
    cout.flush();
  }

  g.clear();

  auto weights = this->m_spliiter->assignWeights(this->m_m, this->m_config);

  if (!m_config.quite)
    cout << " done in " << (float) (clock() - start) / CLOCKS_PER_SEC << " s"
        << endl;
}

void Unfolder::buildMST(GRAPH& g)
{
  auto start = clock();

  priority_queue<DualGraphEdge, vector<DualGraphEdge>, greater<DualGraphEdge>> q;

  if (!m_config.quite) {
    cout << "- building MST...";
    cout.flush();
  }

  this->m_parents.clear();
  //this->m_selected_edges.clear();
  //this->m_selected_edges=vector<bool>(m_m->e_size,false);
  memset(this->m_selected_edges, false, sizeof(bool)*m_m->e_size);
  this->m_fold_edges.clear();
  this->m_ordered_face_list.clear();
  this->m_ordered_crease_list.clear();

  m_parents.resize(m_m->t_size);

  //auto in_the_tree = set<uint>();
  DisjointSets DS(m_m->t_size); //Disjoint set

  // base_face is in the tree
  m_ordered_face_list.push_back(m_base_face_id);
  m_parents[m_base_face_id] = -1;

  // insert all edges start from base face to the queue
  for (const auto& de : g[m_base_face_id]) {
    q.push(DualGraphEdge(m_base_face_id, de.first, de.second));
  }

  float total_weight = 0;

  // total t_size - 1 (- boundary_edges) edges in the tree
  for (uint k = 1; k < m_m->t_size; k++)
  {
    //auto e = q.top();
    DualGraphEdge e;

    // find the maximum weight edge with one node in the tree, but another not
    while (q.empty() == false) {
      e = q.top();
      q.pop();

      // both node are already in the tree
      // not valid, remove that edge
      if(DS.find(e.fid1)==DS.find(e.fid2)) continue;

      break;
    }

    if (e.fid1 < 0 || e.fid2 < 0) {
      break;
    }

    // add directed edge to selected edge set
    //this->m_selected_edges.insert(make_pair(e.fid1, e.fid2));

    {
      int eid = this->m_m->getEdgeIdByFids(e.fid1, e.fid2);
      assert(eid>=0 && eid<this->m_m->e_size);
      this->m_ordered_crease_list.push_back(eid);
      this->m_fold_edges.insert(eid);
      this->m_selected_edges[eid]=true;
    }

    total_weight += e.weight;

    int tree_set=DS.find(m_base_face_id);
    auto s = (DS.find(e.fid1)==tree_set) ? e.fid1 : e.fid2;
    auto t = s == e.fid1 ? e.fid2 : e.fid1;
    tree_set = DS.unite(DS.find(s),DS.find(t));
    this->m_parents[t] = s;
    this->m_ordered_face_list.push_back(t);

    // add new edges
    for (const auto& de : g[t])
    {
      // if the other face is already in the tree
      if (DS.find(de.first)==tree_set) continue;
      q.push(DualGraphEdge(t, de.first, de.second));
    }
  }//end for

  if (!m_config.quite) {
    cout << " total weight = " << total_weight << endl;
    cout << " fold edges = " << m_fold_edges.size() << endl;
    cout << " done in " << (float) (clock() - start) / CLOCKS_PER_SEC << " s"   << endl;
  }
}

double Unfolder::rebuildTree(int base_face, const bool * selected_edges) {

  //reset
  this->m_parents.clear();
  this->m_ordered_face_list.clear();
  this->m_parents.resize(m_m->t_size);

  // auto in_the_tree = set<uint>();
  // in_the_tree.insert(base_face);

  //init
  this->m_ordered_face_list.push_back(base_face);
  this->m_parents[base_face] = -1;
  this->m_m->tris[base_face].path_len = 0;
  this->m_max_path_len = -1;
  this->m_avg_path_len = -1;

  list<int> open;
  open.push_back(base_face);
  auto sum_path_len = 0.0;

  //flooding
  while(!open.empty())
  {
    int fid=open.front();
    open.pop_front();
    triangle & tri=this->m_m->tris[fid];

    for(short i=0;i<3;i++)
    {
      int eid=tri.e[i];
      if(!m_selected_edges[eid]) continue; //not a crease
      int ofid=this->m_m->edges[eid].otherf(fid);
      if(this->m_parents[fid]==ofid) continue; //ofid is parent

      //update info of ofid
      this->m_parents[ofid] = fid;
      this->m_ordered_face_list.push_back(ofid);
      auto path_len=this->m_m->tris[ofid].path_len = this->m_m->tris[fid].path_len + 1;
      this->m_max_path_len = max(this->m_max_path_len, path_len);
      sum_path_len += path_len;

      //add to open to propagate
      open.push_back(ofid);
    }//end for i

  }//end while

  //TODO: JML, this is brute force....should be improved
  // while (selected_edges.size()) {
  //   set<pair<int, int>> to_remove;
  //   for (const auto& e : selected_edges) {
  //     if (!in_the_tree.count(e.first) && !in_the_tree.count(e.second))
  //       continue;
  //
  //     int pf = in_the_tree.count(e.first) ? e.first : e.second;
  //     int f = pf == e.first ? e.second : e.first;
  //
  //     in_the_tree.insert(f);
  //     this->m_parents[f] = pf;
  //     this->m_ordered_face_list.push_back(f);
  //
  //     this->m_m->tris[f].path_len = this->m_m->tris[pf].path_len + 1;
  //
  //     this->m_max_path_len = max(this->m_max_path_len,
  //         this->m_m->tris[f].path_len);
  //
  //     sum_path_len += this->m_m->tris[f].path_len;
  //
  //     to_remove.insert(e);
  //   }
  //
  //   for (const auto& e : to_remove) {
  //     selected_edges.erase(e);
  //   }
  // }

  assert(this->m_ordered_face_list.size() == this->m_m->t_size);

  this->m_avg_path_len = sum_path_len / this->m_m->t_size;

  cout << " - Base face = " << base_face << " avg_path_length = "
       << m_avg_path_len << "\r" << flush;

  return this->m_avg_path_len;
}

void Unfolder::findBestBaseFace() {

  cerr << "- Finding best base face... " << endl;

  int best_base_face = -1;
  double min_avg_path_len = FLT_MAX;

  for (int i = 0; i < m_m->t_size; ++i)
  {
    this->rebuildTree(i, this->m_selected_edges);
    if (m_avg_path_len < min_avg_path_len) {
      min_avg_path_len = m_avg_path_len;
      best_base_face = i;
    }
  }

  this->m_config.baseface = this->m_base_face_id = best_base_face;

  this->rebuildTree(best_base_face, this->m_selected_edges);

  this->rebuildModel();

  cerr << "- Best base face = " << best_base_face << " avg_path_length = "
       << min_avg_path_len << endl;

}

void Unfolder::initUnfold() {
  auto start = clock();

  if (!m_config.quite) {
    cout << "- Initializing unfolding" << endl;
    cout << "  - vertices = " << m_m->v_size		// number of vertices
        << " edges = " << (m_m->e_size - m_m->e_boundary_size)// number of non-bounday edges / total edges
        << "/" << m_m->e_size << " faces = " << m_m->t_size << endl;// number of total facet
  }

  // use random base face
  if (this->m_config.random_baseface)
    m_base_face_id = (int) (mathtool::drand48() * m_m->t_size);
  else if (this->m_config.baseface >= 0) {
    if (this->m_config.baseface >= this->m_m->t_size) {
      std::cerr << "! Warning: base face out of range! baseface="<<this->m_config.baseface<<", t size="<<this->m_m->t_size<< endl;
      m_base_face_id = 0;
      //exit(-1);
    }
    m_base_face_id = this->m_config.baseface;
  } else {
    m_base_face_id = 0;
  }

  if (!m_config.quite)
    cout << "  - base face = " << m_base_face_id << endl;

  this->m_org.clear();

  // copy vertex coordinates
  for (auto i = 0; i < this->m_m->t_size; i++) {
    const auto& f = this->m_m->tris[i];
    const auto& vs = this->m_m->vertices;
    auto v0 = std::make_pair(f.v[0], Vector3d(vs[f.v[0]].p.get()));
    auto v1 = std::make_pair(f.v[1], Vector3d(vs[f.v[1]].p.get()));
    auto v2 = std::make_pair(f.v[2], Vector3d(vs[f.v[2]].p.get()));

    m_org.push_back( { v0, v1, v2 });
  }

// copy to unfolded
  m_unfolded = m_org;

// reset flat edges
  m_flat_edges = 0;

// shared edges
  m_shared_edge.clear();
  for (auto i = 0; i < m_m->e_size; i++) {
    const auto& edge = m_m->edges[i];
    const auto fid1 = edge.fid[0];
    const auto fid2 = edge.fid[1];
    m_shared_edge[fid1][fid2] = m_shared_edge[fid2][fid1] = i;

    if (edge.folding_angle == 0)
      m_flat_edges++;
  }

  if (!m_config.quite)
    cout << "  - Done in " << (float) (clock() - start) / CLOCKS_PER_SEC << " s"
        << endl;
}

void Unfolder::alignModel() {
// align based_face's normal with Y axis
  auto n0 = this->m_m->tris[m_base_face_id].n.normalize();
  // work around
  auto y_axis = Vector3d(0, 1, 0).normalize();
  this->m_rotation_angle = acos(n0 * y_axis);
  this->m_rotation_axis = n0 % y_axis;

// no rotation required
  if (this->m_rotation_angle == 0) {
    this->m_rotation_axis = Vector3d(0, 1, 0);
  } else if (fabs((fabs(this->m_rotation_angle) - PI)) < 1e-6) {
    this->m_rotation_axis = Vector3d(1, 0, 0);
    this->m_rotation_angle = PI;
  }

  auto r_matrix = Matrix4x4::getRotationMatrix(Vector3d(0, 0, 0),
      this->m_rotation_axis, this->m_rotation_angle);

  for (auto& f : this->m_org) {
    for (auto k = 0; k < 3; k++) {
      f[k].second = (r_matrix * f[k].second.from3dto4d()).from4dto3d();
    }
  }

  const auto v0 = this->m_org[m_base_face_id][0].second;

  this->m_transilation = -v0;

  for (auto& f : this->m_org) {
    for (auto k = 0; k < 3; k++) {
      f[k].second = f[k].second + m_transilation;
    }
  }

  this->m_unfolded = this->m_org;
}

void Unfolder::computeUnfolding() {
  if (!m_config.quite) {
    cout << "- Computing unfolding...";
    cout.flush();
  }

  auto start = clock();

  this->linearUnfoldTo(1.0);

  if (!m_config.quite)
    cout << " Done in " << (float) (clock() - start) / CLOCKS_PER_SEC << " s"
        << endl;
}

const Vector3d& Unfolder::getUnfoldedVertex(uint fid, uint org_vid) const {
  const auto& f = this->m_unfolded[fid];
  for (const auto& pair : f) {
    if (pair.first == org_vid)
      return pair.second;
  }

  assert(false);

// remove warning...
  return this->m_unfolded[0][0].second;
}

bool Unfolder::isEdgeCCW(uint fid, uint vid1, uint vid2) const {
  for (auto k = 0; k < 3; k++) {
    const auto fvid1 = m_m->tris[fid].v[k];
    const auto fvid2 = m_m->tris[fid].v[(k + 1) % 3];

    if (fvid1 == vid1 && fvid2 == vid2)
      return true;
  }

  return false;
}

// find the boundary of the net
void Unfolder::findBoundary() {
// clear the result
  this->m_boundary.clear();

  unordered_set<int> attached_vertices;

// an edge from vid1 - vid2
  int vid1 = INT_MAX;
  int vid2 = INT_MAX;

  for (int fid = 0; fid < this->m_m->t_size; ++fid) {
    for (int i = 0; i < 3; ++i) {
      const auto s = this->m_fs[fid][i];
      const auto t = this->m_fs[fid][(i + 1) % 3];

      const auto edge = make_pair(s, t);

      if (!this->m_crease_lines.count(edge)) {
        // we find the first edge
        vid1 = s;
        vid2 = t;
        break;
      }
    }

    // found an boundary edge
    if (vid1 != INT_MAX && vid2 != INT_MAX) {
      cout << " - initial boundary edge = " << vid1 << " -> " << vid2
          << " at face = " << fid << endl;
      break;
    }
  }

  m_boundary.push_back(vid1);
  m_boundary.push_back(vid2);

  attached_vertices.insert(vid1);
  attached_vertices.insert(vid2);

// we haven't find a loop
  while (vid2 != vid1) {
    bool found = false;
    // loop all adjacent faces of vid2
    for (auto fid : m_vfmap[vid2]) {
      for (int i = 0; i < 3; ++i) {
        auto s = this->m_fs[fid][i];
        auto t = this->m_fs[fid][(i + 1) % 3];

        // next edge must started with vid2
        if (s != vid2)
          continue;

        auto edge = make_pair(s, t);

        // it's the boundary
        if (!this->m_crease_lines.count(edge)) {
          found = true;

          // we find the edge
          m_boundary.push_back(t);
          attached_vertices.insert(t);

          vid2 = t;

          break;
        }
      }

      if (found)
        break;
    }
  }

  if (!m_config.quite) {
    cout << "boundary:" << endl;
    for (auto vid : m_boundary)
      cout << vid << " ";
    cout << endl;
  }
}

// find single path for all the creases
void Unfolder::findSinglePath() {
  this->m_single_path.clear();

  int vid1, vid2;

  bool found = false;

// find the start vertex and start edge
  for (auto i = 0; !found && i < this->m_vs.size(); ++i) {
    // loop all adjacent faces
    for (auto fid : this->m_vfmap[i]) {
      for (auto k = 0; k < 3; ++k) {
        auto s = this->m_fs[fid][k];
        auto t = this->m_fs[fid][(k + 1) % 3];

        // an edge must start with vid1
        if (s != i)
          continue;

        auto edge = make_pair(s, t);

        if (this->m_crease_lines.count(edge)) {
          vid1 = s;
          vid2 = t;
          found = true;
          break;
        }
      } // end for k
    } // end for fid
  } // end for i

  unordered_set<int> vistied_vs;
  set<pair<int, int>> visited_edges;

  this->findSinglePath(vid1, vistied_vs, visited_edges, this->m_single_path);

  cout << "single path = " << endl;
  for (auto vid : m_single_path)
    cout << vid << " ";
  cout << endl;
}

// travel from vid
void Unfolder::findSinglePath(int vid, unordered_set<int>& visited_vids,
    set<pair<int, int>>& visited_edges, vector<int>& path) {
// append current node to the path
  path.push_back(vid);

// mark vid as visited
  visited_vids.insert(vid);

// loop adjacent faces
  for (const auto fid : this->m_vfmap[vid]) {
    // loop all edges
    for (auto k = 0; k < 3; ++k) {
      const auto s = this->m_fs[fid][k];
      const auto t = this->m_fs[fid][(k + 1) % 3];

      // an edge must start with vid
      if (s != vid)
        continue;

      auto edge = make_pair(s, t);

      // is a crease line
      if (this->m_crease_lines.count(edge)) {
        // direct goto if edge not visited but vertex visited
        if (visited_vids.count(t) && !visited_edges.count(edge)) {
          visited_edges.insert(edge);
          path.push_back(t);
        }
        // recursive goto the vid
        else {
          visited_edges.insert(edge);
          this->findSinglePath(t, visited_vids, visited_edges, path);
        }

        // move back
        path.push_back(-vid - 1);
      }
    } // end for k
  } // end for fid
}

bool Unfolder::getIntersectionOfPlane(const int fid1, const Vector3d& v1,
    const Vector3d& v2, Vector3d& outPoint) {
  const double EP = 1e-3;

// the normal of the edge
  Vector3d vl = (v1 - v2).normalize();

// get 3 vertices of the face
  Vector3d p0 = this->m_unfolded[fid1][0].second;
  Vector3d p1 = this->m_unfolded[fid1][1].second;
  Vector3d p2 = this->m_unfolded[fid1][2].second;

// two vector on the face
  Vector3d vf0 = (p1 - p0);
  Vector3d vf1 = (p2 - p0);

// connect the end points of the edge to one point on the plane
  Vector3d vl1 = (v1 - p0);
  Vector3d vl2 = (v2 - p0);

// normal vector of the face
  Vector3d fno = (vf0 % vf1).normalize();

  double dp1 = vl1 * fno;
  double dp2 = vl2 * fno;

// line segment are on the same side from the plane
  if (sign(dp1) == sign(dp2) || (fabs(dp1) < EP && fabs(dp2) < EP))
    return false;

  double dot_prodcut = (vl * fno);

  if (fabs(dot_prodcut) < EP) {
    //the line is parallel to the plane
    //???
    return false;
  } else {
    double t = -vl1 * fno / dot_prodcut;
    if (fabs(t) < EP)
      return false;

    //calculate the intersection point
    outPoint[0] = v1[0] + vl[0] * t;
    outPoint[1] = v1[1] + vl[1] * t;
    outPoint[2] = v1[2] + vl[2] * t;

    Vector3d pOut = Vector3d(outPoint[0], outPoint[1], outPoint[2]);

    for (uint i = 0; i < 3; i++) {
      Vector3d p = this->m_unfolded[fid1][i].second;
      if ((p - pOut).norm() < EP) {
        return false;
      }
    }
    return true;
  }
}

bool Unfolder::pointInTriangle(const int fid1, const Vector3d& p) {
  Vector3d p0 = this->m_unfolded[fid1][0].second;
  Vector3d p1 = this->m_unfolded[fid1][1].second;
  Vector3d p2 = this->m_unfolded[fid1][2].second;

  Vector3d v0 = (p1 - p0);
  Vector3d v1 = (p2 - p0);
  Vector3d v2 = (p - p0);

  double dot00 = v0 * v0;
  double dot01 = v0 * v1;
  double dot02 = v0 * v2;
  double dot11 = v1 * v1;
  double dot12 = v1 * v2;

  double invDenom = 1 / (dot00 * dot11 - dot01 * dot01);
  double u = (dot11 * dot02 - dot01 * dot12) * invDenom;
  double v = (dot00 * dot12 - dot01 * dot02) * invDenom;

// Check if point is in triangle
  if ((u >= 0) && (v >= 0) && (u + v < 1)) {
    return true;
  }

  return false;
}

bool Unfolder::hasIntersection(const int fid1, const int fid2) {
  bool intersection = false;

  const auto face1 = this->m_m->tris[fid1];
  const auto face2 = this->m_m->tris[fid2];

  Vector3d p;

  for (uint i = 1; i <= 3; i++) {
    const Vector3d& v1 = this->m_unfolded[fid2][i - 1].second;
    const Vector3d& v2 = this->m_unfolded[fid2][i % 3].second;

    // get the intersection of the line contains the edge and the plane contains the face
    if (this->getIntersectionOfPlane(fid1, v1, v2, p)) {
      // check whether the intersection is in the face
      if (this->pointInTriangle(fid1, p)) {
        intersection = true;
      }
    }
  }

  return intersection;
}

uint Unfolder::checkOverlap_itree() {
  // move this to other file to speed up the compilation
  ITreeChecker checker(this);
  return checker.checkOverlapping(this->m_unfolded, this->m_config);
}

bool Unfolder::checkOverlap(uint i, uint j)
{
  bool intersected = false;
  double pp[2] = { 0.0, 0.0 };

  // two faces shared an edge in the unfolding (not in the original mesh).
  int eid=this->m_m->getEdgeIdByFids(i,j);
  if(eid>=0 && eid<this->m_m->e_size) if (this->m_selected_edges[eid]) return false;

  for (int p = 0; p < 3 && (!intersected); p++) {

    const auto& vi1 = m_unfolded[i][p].second;
    const auto& vi2 = m_unfolded[i][(p + 1) % 3].second;

    const double a[2] = { vi1[0], vi1[2] };
    const double b[2] = { vi2[0], vi2[2] };

    for (int q = 0; q < 3 && (!intersected); q++) {

      const auto& vj1 = m_unfolded[j][q].second;
      const auto& vj2 = m_unfolded[j][(q + 1) % 3].second;

      if ((vi1 == vj1 && vi2 == vj2) || (vi1 == vj2 && vi2 == vj1))
        continue;

      const double c[2] = { vj1[0], vj1[2] };
      const double d[2] = { vj2[0], vj2[2] };

      char r = SegSegInt<double>(a, b, c, d, pp);

      if (r == '1') {
        intersected = true;
      }
    }
  }

  if (intersected) {
    this->m_m->tris[i].overlapped = true;
    this->m_m->tris[j].overlapped = true;
    if (this->m_config.record_overlap) {
      this->m_overlap_pairs[i].insert(j);
      this->m_overlap_pairs[j].insert(i);
    }
  }

  return intersected;
}

bool Unfolder::checkOverlapNew(uint i, uint j)
{
  bool intersected = false;
  // two faces shared an edge in the unfolding (not in the original mesh).
  int eid=this->m_m->getEdgeIdByFids(i,j);
  if(eid>=0 && eid<this->m_m->e_size) if (this->m_selected_edges[eid]) return false;

  // if (this->m_selected_edges.count(make_pair(i, j))
  //     || this->m_selected_edges.count(make_pair(j, i)))
  //   return intersected;

  double p1[2] = { m_unfolded[i][0].second[0], m_unfolded[i][0].second[2] };
  double q1[2] = { m_unfolded[i][1].second[0], m_unfolded[i][1].second[2] };
  double r1[2] = { m_unfolded[i][2].second[0], m_unfolded[i][2].second[2] };

  double p2[2] = { m_unfolded[j][0].second[0], m_unfolded[j][0].second[2] };
  double q2[2] = { m_unfolded[j][1].second[0], m_unfolded[j][1].second[2] };
  double r2[2] = { m_unfolded[j][2].second[0], m_unfolded[j][2].second[2] };

  intersected = tri_tri_overlap_test_2d(p1, q1, r1, p2, q2, r2);
  if (intersected) {
    this->m_m->tris[i].overlapped = true;
    this->m_m->tris[j].overlapped = true;
    if (this->m_config.record_overlap) {
      this->m_overlap_pairs[i].insert(j);
      this->m_overlap_pairs[j].insert(i);
    }
  }

  return intersected;
}

//TODO: this is too slow
int Unfolder::checkOverlaps() {

  ++this->m_check_overlapping_calls;

  const int F = this->m_unfolded.size();

  if (!m_config.quite) {
    cout << "- checking overlapping... ";
    cout.flush();
  }

// backup the current status
  MESH unfoled_back { this->m_unfolded };

// shrink the triangles only for checking overlapping...
  if (this->m_config.shrink) {
    this->shrink();
  }

  auto start = clock();

  uint count = 0;

  if (this->m_config.record_overlap) {
    this->m_overlap_pairs.clear();
    m_overlap_pairs.resize(F);
  }

  for (int i = 0; i < F; i++)
    this->m_m->tris[i].overlapped = false;

  if (!m_config.use_rapid) {

    bool checked = false;
    if (F > 300) //if the number of faces is large, use more advanced data structure
    {
      count = checkOverlap_itree();
      if (count != UINT_MAX)
        checked = true; //good
    }

    if (checked == false) {

      // Build AABB for the unfolding
      // vector<Box2d> boxes(F);
      //
      // for (int i = 0; i < F; ++i) {
      //   Vector2d p[3];
      //   for (short j = 0; j < 3; ++j)
      //     p[j] = {m_unfolded[i][j].second[0], m_unfolded[i][j].second[2]};
      //   boxes[i].setFromPoints(p);
      // }

      for (int i = 0; i < F; i++) {
        for (int j = i + 1; j < F; j++) {

       //   if (!boxes[i].intersect(boxes[j]))
        //    continue;

          if (checkOverlapNew(i, j)) {
          //if (checkOverlap(i, j)) {
            count++;
          }
        } //end for j
      } //end for i
    } //end if (checked == false)

  } else {
    //use rapid
    count = this->m_cd->hasCollision();
  }

  const auto total_time = (clock() - start) * 1.0 / CLOCKS_PER_SEC;

  this->m_total_check_overlapping_time += (clock() - start);

  if (!m_config.quite)
    cout << " count = " << count << " Done in " << total_time << " s" << endl;

// recovery
  this->m_unfolded.swap(unfoled_back);

  this->m_is_flattened = (count == 0);

  m_last_overlap_count = count;

  return count;
}

// count local overlaps due to insufficient count for hyperbolic vertex
int Unfolder::checkLocalOverlaps()
{
  int overlaps = 0;
  //int concave_vertex = 0;

  for (auto i = 0; i < this->m_m->v_size; ++i) {
    const auto& v = this->m_m->vertices[i];

    // only consider hyperbolic vertices
    if (!v.hyperbolic)
      continue;

    //concave_vertex++;

    //int cuts = 0;
    for (const auto eid : v.m_e) {
      const edge e = this->m_m->edges[eid];
      if (e.fid.size() == 1)
        continue;
      const int fid1 = e.fid[0];
      const int fid2 = e.fid[1];

      if (this->isCutEdge(eid) && this->m_overlap_pairs[fid1].count(fid2)) {
        overlaps++;
      }
    }//end for eid
  }//end for i

  return overlaps;
}

// collision is detected based on penetration
bool Unfolder::checkCollision() {
  const int F = this->m_unfolded.size();

  if (!m_config.quite) {
    cout << "- checking collision... ";
    cout.flush();
  }

  auto start = clock();

  int count = 0;

  for (int i = 0; i < F; i++)
    this->m_m->tris[i].overlapped = false;

  for (int i = 0; i < F; i++) {
    for (int j = i+1; j < F; j++) {
      // if (i == j)
      //   continue;

      // two faces shared an edge in the unfolding (not in the original mesh).
      int eid=this->m_m->getEdgeIdByFids(i,j);
      if(eid>=0 && eid<this->m_m->e_size) if (this->m_selected_edges[eid]) continue;

      bool intersected = this->hasIntersection(i, j);

      if (intersected) {
        ++count;
        this->m_m->tris[i].overlapped = true;
        this->m_m->tris[j].overlapped = true;
      }
    }//end for j
  }//end for i

  if (!m_config.quite)
    cout << " count = " << count << " done in "
        << (float) (clock() - start) / CLOCKS_PER_SEC << " s" << endl;

  return count > 0;
}

void Unfolder::rebuildModel() {

  if (m_config.no_rebuild) {
    cerr << "! Warning! NO rebuild, skipped! [Unfolder::rebuildModel]" << endl;
    return;
  }

  if (!m_config.quite)
    cerr << "- Rebuilding model..." << endl;

  this->initUnfold();

  this->alignModel();

  // Step 0: unfold this to flat
  this->unfoldTo(1.0);

  // Step 1: copy the model...
  this->m_net = *this->m_m;

  // Step 2: cut the model.

  // Step 2.1: collect all dual edges  {fid1, fid2}.
  set<pair<int, int>> dual_edges;
  uint eid=-1;
  for (const edge& e : m_m->edges) {
    eid++;
    if (e.type == 'b') continue;
    if (m_selected_edges[eid]) continue; //crease edge, ignore
    dual_edges.insert(make_pair((int) e.fid[0], (int) e.fid[1]));
  }

  // Step 2.2: remove all fold dual edges.
  // for (const auto& p : m_selected_edges) {
  //   dual_edges.erase(make_pair(p.first, p.second));
  //   dual_edges.erase(make_pair(p.second, p.first));
  // }

  // cut edge {fid1, fid2}
  // ?? not sure why this is needed
  vector<pair<int, int>> cut_dual_edges(dual_edges.begin(), dual_edges.end());

  // Step 2.3: cut the mesh
  mesh_cutter::CutMesh(&this->m_net, cut_dual_edges);

  // Step 3: update the coordinates for each face in the new model
  this->sync(&this->m_net);

//  cerr << " done in " << endl;
}

/// Sync correct coordinates with an unfolding.
/// Will modify model* unfolding.
void Unfolder::sync(model* unfolding) const {
  for (int i = 0; i < m_m->t_size; ++i) {
    const triangle& f = unfolding->tris[i];
    for (int j = 0; j < 3; ++j) {
      unfolding->vertices[f.v[j]].p.set(this->m_unfolded[i][j].second.get());
    }
  }
}

// add tabs on the net
//JML: deprecated. this is implemented now in svg writer
#if 0
void Unfolder::addTabs() {
// add a tab for each cut edge in the original model
// tab width = edge width
// tab height = edge width * 15%
// tab angle = try [45, 30, 15]

// steps:
// 1. find cut/boundary edges to add tab
// 2. add two triangle for each edge

  int count_added = 0;
  int count_total = 0;

  for (const auto& ee : m_eemap) {
    set<int> vids;
    for (const auto pair : ee.second) {
      vids.insert(pair.first);
      vids.insert(pair.second);
    }

    bool boundary_edge = this->m_m->edges[ee.first].fid.size() == 1;

    if (boundary_edge) {
      cout << "boundary edge = " << ee.first << endl;
    }

    // a cut edge has at least 3 vertices in the net
    if (vids.size() >= 3 || boundary_edge) {
      // find an cut edge
      ++count_total;

      if (this->addTabForCutEdge(ee.first))
      ++count_added;
    }
  }

  cout << "\n\n";
  cout << "Tabs added = " << count_added << "/" << count_total << endl;
}
#endif

int Unfolder::getNetFidByVids(const int vid1, const int vid2) {
  set<int> vids = { vid1, vid2 };

  for (int i = 0; i < this->m_fs.size(); ++i) {
    const auto face = this->m_fs[i];
    int count = 0;

    for (int j = 0; j < 3; ++j) {
      if (vids.count(face[j]))
        ++count;
    }

    if (count == 2)
      return i;
  }

  return -1;
}

//JML: deprecated. this is implemented now in svg writer
#if 0
/// add one tab for each cut edge with original edge id
bool Unfolder::addTabForCutEdge(const int org_eid) {
  for (const auto& e : this->m_eemap[org_eid]) {
    const int vid1 = e.first;
    const int vid2 = e.second;
    const int fid = getNetFidByVids(vid1, vid2);

    if (this->addTabForNewEdge(vid1, vid2, fid)) {
      return true;
    }
  }

  cerr << "can't add tab on to edge " << org_eid << endl;

  return false;
}

/// add one tab on an edge of the net
bool Unfolder::addTabForNewEdge(const int vid1, const int vid2, const int fid) {
// see: https://www.dropbox.com/s/16yrtfxnkomubwi/mesh_unfolder_tab.png?dl=0

// get the coordinates of the edge for adding tab
  const auto& p1 = this->m_vs[vid1];
  const auto& p2 = this->m_vs[vid2];
// vector of the edge
  auto edge = p2 - p1;
// length of the edge
  const auto edge_len = edge.norm();
// normalize the edge
  edge = edge.normalize();

//TODO, should also be based on the size of the entire net
  const auto min_net_size = min(m_net_size[0], m_net_size[2]);
  const auto tab_heights = {min_net_size * 0.05, min_net_size * 0.04,
    min_net_size * 0.03};

  const vector<double> angles = {45.0, 30.0, 20.0, 10.0};

// try different combinations...
  for (const auto tab_height : tab_heights) {
    for (const auto angle3 : angles) {
      for (const auto angle4 : angles) {
        const auto m3 = Matrix4x4::getRotationMatrixY(DegToRad(-angle3));
        const auto m4 = Matrix4x4::getRotationMatrixY(DegToRad(angle4));
        const auto len3 = tab_height / sin(DegToRad(angle3));
        const auto len4 = tab_height / sin(DegToRad(angle4));

        // |p3p4|
        const auto len34 = edge_len
        - (tab_height / tan(DegToRad(angle3))
            + tab_height / tan(DegToRad(angle4)));

        if (len34 < 0 || len34 * 2 < edge_len)
        continue;

        // compute the vertices of the tab
        const auto p3 = (m3 * edge.from3dto4d()).from4dto3d() * len3 + p1;
        const auto p4 = (m4 * (-edge).from3dto4d()).from4dto3d() * len4 + p2;

        //check collision, latter
        if (this->isOverlapWithNet(p1, p3, p4, fid)
            || this->isOverlapWithNet(p1, p4, p2, fid)) {
          continue;
        }

        // add tab to the net...
        // 1. add vertices
        // 2. add creases
        // 3. add faces
        // 3. update ordered face list

        // create new vertices
        const int vid3 = this->m_vs.size();
        const int vid4 = this->m_vs.size() + 1;

        this->m_vs.push_back(p3);
        this->m_vs.push_back(p4);

        // create new faces
        int fid134 = this->m_fs.size();
        int fid142 = this->m_fs.size() + 1;

        // face 1 3 4
        this->m_fs.insert(this->m_fs.end(), {vid1, vid3, vid4});

        // face 1 4 2
        this->m_fs.insert(this->m_fs.end(), {vid1, vid4, vid2});

        // create new crease lines
        int cid12 = this->m_cs.size();
        int cid14 = this->m_cs.size() + 1;

        auto c12 = Crease(vid1, vid2, fid142, fid, PI);
        auto c14 = Crease(vid1, vid4, fid134, fid142, 0);
        this->m_cs.push_back(c12);
        this->m_cs.push_back(c14);

        this->m_crease_lines[make_pair(vid1, vid2)] = cid12;
        this->m_crease_lines[make_pair(vid1, vid4)] = cid14;

        this->m_ordered_face_list.push_back(fid142);
        this->m_ordered_face_list.push_back(fid134);

        return true;

      }   // for angle3
    }   // for angle4
  }   // for heights

  return false;
}
#endif

/// check whether a triangle overlap with existing net (included attached tabs)
bool Unfolder::isOverlapWithNet(const Vector3d& p1, const Vector3d& p2,
    const Vector3d& p3, const int fid) {
// output, intersection location
  double pp[2] = { 0.0, 0.0 };

  const auto tab = vector<Vector3d> { p1, p2, p3 };

  int cur_fid = 0;
  for (const auto& f : this->m_fs) {
    // do not test against parent face...
    if (cur_fid++ == fid)
      continue;

    for (int i = 1; i <= 3; ++i) {
      const auto& vi1 = this->m_vs[f[i - 1]];
      const auto& vi2 = this->m_vs[f[i % 3]];

      for (int j = 1; j <= 3; ++j) {
        const auto& vj1 = tab[j - 1];
        const auto& vj2 = tab[j % 3];

        if ((vi1 == vj1 && vi2 == vj2) || (vi1 == vj2 && vi2 == vj1))
          continue;

        const double a[2] = { vi1[0], vi1[2] };
        const double b[2] = { vi2[0], vi2[2] };

        const double c[2] = { vj1[0], vj1[2] };
        const double d[2] = { vj2[0], vj2[2] };

        char r = SegSegInt<double>(a, b, c, d, pp);

        if (r == '1') {
          cout << "has overlap!!!!" << endl;
          return true;
        }
      }
    }
  }

  return false;
}

bool Unfolder::isCutEdge(int fid1, int fid2) const {
  uint eid = this->m_m->getEdgeIdByFids(fid1, fid2);
  return isCutEdge(eid);
  //return !this->m_selected_edges[eid];
  //return !this->m_selected_edges.count( { fid1, fid2 })
  //    && !this->m_selected_edges.count( { fid2, fid1 });
}

bool Unfolder::isCutEdge(int eid) const {

  const edge & e = this->m_m->edges[eid];

  // A border edge can not be a cut edge.
  if (e.type == 'b')
    return false;

  return !this->m_selected_edges[eid];
}

vector<float> Unfolder::optimizeUnfolding() {

  if (!m_config.quite)
    cerr << "- Optimizing unfolding..." << endl;

  this->m_optimization_seqs.clear();

  // {(length eid)}
  vector<pair<float, int>> edge_lengths(m_m->e_size);
  vector<float> el(m_m->e_size);

  // Use current ones as the best ones
  vector<float> best_weights = this->m_weights;
  vector<float> my_weights = this->m_weights;
  float best_cut_length = this->getTotalCutLength();
  set<uint> best_fold_edges = this->m_fold_edges;

  if (!m_config.quite)
    cerr << "- init best cut_length = " << best_cut_length << endl;

  // Compute edge lengths
  for (int i = 0; i < m_m->e_size; ++i) {

    const edge& e = this->m_m->edges[i];
    const Point3d& v0 = this->m_m->vertices[e.vid[0]].p;
    const Point3d& v1 = this->m_m->vertices[e.vid[1]].p;

    edge_lengths[i].first = el[i] = (float) (v1 - v0).norm();
    edge_lengths[i].second = i;

  }

  const float INF_WEIGHT = 1e6;

#if 0
  for (int i = 0; i < m_m->t_size; ++i) {
    const triangle& f = this->m_m->tris[i];

    // { edge_length, eid }
    vector<pair<float, int>> es;

    for (int j = 0; j < 3; ++j) {
      es.push_back(std::make_pair(el[f.e[j]], (int) f.e[j]));
    }

    std::sort(es.begin(), es.end());

    vector<float> weights = best_weights;

    if (best_fold_edges.count(es[0].second)
        && (!best_fold_edges.count(es[1].second)
            || !best_fold_edges.count(es[2].second))) {
      weights[es[0].second] = INF_WEIGHT; // cut the shortest edge of the triangle.
      int overlaps = this->buildFromWeights(weights);
      if (overlaps != 0)
      continue;

      float cut_length = this->m_spliiter->getTotalCutLength(this);

      if (cut_length < best_cut_length) {
        best_cut_length = cut_length;
        best_fold_edges = this->m_fold_edges;
        best_weights[es[0].second] = INF_WEIGHT;

        auto hull = this->buildConvexHull2D();
        auto area = ConvexHull2DArea(hull);

        // Add to evolve seq
        this->m_evolving_seqs.push_back(best_weights);
        cerr << "- new cut length = " << best_cut_length << " hull area = "
        << area << endl;
      }
    }

  }
#endif

#if 1
  // Sort by edge length
  std::sort(edge_lengths.begin(), edge_lengths.end(),
      std::greater<pair<float, int>>());

  for (int i = 0; i < edge_lengths.size(); ++i) {
    int eid = edge_lengths[i].second;

    vector<float> weights = my_weights;

    weights[eid] = -INF_WEIGHT + i; // The edge will be kept
    int overlaps = this->buildFromWeights(weights);
    if (overlaps != 0)
      continue;

    my_weights[eid] = weights[eid];

    float cut_length = this->getTotalCutLength();

    if (cut_length < best_cut_length) {
      best_cut_length = cut_length;
      best_fold_edges = this->m_fold_edges;
      best_weights = my_weights;

      auto hull = this->buildConvexHull2D();
      auto area = ConvexHull2DArea(hull);

      UnfoldingState state(best_weights);
      state.AddProperty("cut length", std::to_string(cut_length));
      state.AddProperty("hull area", std::to_string(area));

      // Add to evolve seq
      this->m_optimization_seqs.push_back(state);
      if (!m_config.quite)
        cerr << "- new_best " << state << endl;
    }
  }
#endif

  return best_weights;
}

float Unfolder::getTotalCutLength() const {
  float total_selected_edge_length = 0.0f;

  for (const auto& eid : this->getFoldEdges())
    total_selected_edge_length += m_m->edges[eid].length;

  return m_m->total_edge_length - total_selected_edge_length;
}

float Unfolder::getHullArea() const {
  auto hull = this->buildConvexHull2D();
  double hull_area = masc::util::ConvexHull2DArea(hull);
  return hull_area;
}

vector<Vector3d> Unfolder::buildConvexHull2D() const {
  vector<Vector3d> input;
  vector<Vector3d> output;

  for (const auto& f : this->m_unfolded) {
    for (int i = 0; i < 3; ++i)
      input.push_back(f[i].second);
  }

  masc::util::ConvexHull2D(input, &output);

  return output;
}

// =================================================================
// public access
// =================================================================

model * Unfolder::getModel() const {
  return this->m_m;
}

const model * Unfolder::getNet() const {
  return &m_net;
}

const Config& Unfolder::getConfig() const {
  return this->m_config;
}

const MESH& Unfolder::getUnfolded() const {
  return this->m_unfolded;
}

const MESH& Unfolder::getOrg() const {
  return this->m_org;
}

const string& Unfolder::getFilename() const {
  return this->m_config.filename;
}

const Vector3d& Unfolder::getRotationAxis() const {
  return this->m_rotation_axis;
}

/// angle in degree
const double Unfolder::getRotationAngle() const {
  return this->m_rotation_angle * 180.0 / PI;
}

const Vector3d& Unfolder::getTranslation() const {
  return this->m_transilation;
}

/// get model color
const Vector3d& Unfolder::getColor() const {
  return this->m_color;
}

// set model color
void Unfolder::setColor(double r, double g, double b) {
  this->m_color.set(r, g, b);
}

const bool Unfolder::isFlattened() const {
  return this->m_is_flattened;
}
