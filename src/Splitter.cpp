/*
 * Splitter.cpp
 *
 *  Created on: Feb 6, 2015
 *      Author: zhonghua
 */

#include "Splitter.h"

#include <cfloat>
#include <unordered_set>
#include <algorithm>
using namespace std;

#include "model.h"
#include "unfolder.h"
#include "util/Statistics.h"

#define INF_WEIGHT 1e6

namespace masc {

Splitter::Splitter() {
  this->m_max_edge_length = FLT_MIN;
  this->m_min_edge_length = FLT_MAX;
  this->m_total_edge_length = 0.0f;
}

void Splitter::measure(model* m) {

  this->m_total_edge_length = 0.0f;

  // add cache for each model
  for (auto i = 0; i < m->e_size; i++) {
    const auto& edge = m->edges[i];
    const auto edge_len = edge.length;

    this->m_max_edge_length = std::max(this->m_max_edge_length, edge_len);
    this->m_min_edge_length = std::min(this->m_min_edge_length, edge_len);

    this->m_total_edge_length += edge_len;
  }
}

Vector3d Splitter::genRandomUnitVector(const Config& config) {
  // generate a random reference vector or use the given vector
  auto c =
      (config.use_user_vector) ?
          config.user_vector :
          Vector3d(mathtool::drand48(), mathtool::drand48(),
              mathtool::drand48()).normalize();
  if (config.use_user_vector == false) {
    for (short i = 0; i < 3; i++)
      if (mathtool::drand48())
        c[i] = -c[i];
  }

  return c;
}

Splitter* Splitter::createSplitter(const string& heruistic) {
  if (heruistic == "MinimumPerimeter")
    return new MinimumPerimeterSplitter();
  if (heruistic == "MaximumPerimeter")
    return new MaximumPerimeterSplitter();
  if (heruistic == "Random")
    return new RandomSplitter();
  if (heruistic == "RandomEdge")
    return new RandomEdgeSplitter();
  if (heruistic == "FlatTree")
    return new FlatTreeSplitter();
  if (heruistic == "UnFlatTree")
    return new UnFlatTreeSplitter();
  if (heruistic == "SteepestEdge")
    return new SteepestEdgeSplitter();
  if (heruistic == "UnSteepestEdge")
    return new UnSteepestEdgeSplitter();

  cerr << "!Error! Unknown splitter type = " << heruistic << endl;

  return nullptr;
}

Splitter* Splitter::createSplitter(CutHeuristic heuristic) {
  switch (heuristic) {
  case CutHeuristic::MINIMUM_PERIMETER:
    return new MinimumPerimeterSplitter();
  case CutHeuristic::MAXIMUM_PERIMETER:
    return new MaximumPerimeterSplitter();
  case CutHeuristic::FLAT_TREE:
    return new FlatTreeSplitter();
  case CutHeuristic::UNFLAT_TREE:
    return new UnFlatTreeSplitter();
  case CutHeuristic::STEEPEST_EDGE:
    return new SteepestEdgeSplitter();
  case CutHeuristic::UNSTEEPEST_EDGE:
    return new UnSteepestEdgeSplitter();
  case CutHeuristic::RANDOM:
    return new RandomSplitter();
  case CutHeuristic::RANDOM_EDGE:
    return new RandomEdgeSplitter();
  case CutHeuristic::BRUTE_FORCE:
    return new BruteForceSplitter();
  default:
    assert(false);
    break;
  }

  return nullptr;
}

vector<float> Splitter::assignWeights(model* m, const Config& config)
{
  vector<float> weights(m->e_size);

  this->assignWeightsImpl(m, weights, config);

  if(config.less_cuts) {
    for (uint i = 0; i < m->e_size; i++) {
      const edge e = m->edges[i];
      if(e.type == 'd') {
        weights[i] = EdgeWeight::DiagnalEdge;
      } else if(e.type == 'p') {
        weights[i] = EdgeWeight::FlatEdge;
      }
    }
  }

  return weights;
}

//////////////////////////////////////////////////////////////////////

void MinimumPerimeterSplitter::assignWeightsImpl(model *m,
    vector<float>& weights, const Config& config) {

  for (auto i = 0; i < m->e_size; i++) {
    const auto& edge = m->edges[i];
    const auto edge_len = m->edges[i].length;

    float weight = 1.0
        - (edge_len - this->m_min_edge_length)
            / (this->m_max_edge_length - this->m_min_edge_length);

    weights[i] = weight;
  }

}

//////////////////////////////////////////////////////////////////////

void MaximumPerimeterSplitter::assignWeightsImpl(model* m,
    vector<float>& weights, const Config& config) {
  MinimumPerimeterSplitter::assignWeightsImpl(m, weights, config);

  for (auto& w : weights) {
    if (w <= 0)
      continue;
    w = 1.0 - w;
  }
}

//////////////////////////////////////////////////////////////////////

void FlatTreeSplitter::assignWeightsImpl(model* m, vector<float>& weights,
    const Config& config) {

  // generate a random reference vector or use the given vector
  const auto c = this->genRandomUnitVector(config);

  for (auto i = 0; i < m->e_size; i++) {
    const auto& edge = m->edges[i];
    const auto fid1 = edge.fid[0];
    const auto fid2 = edge.fid[1];

    const auto& v1 = m->vertices[edge.vid[0]];
    const auto& v2 = m->vertices[edge.vid[1]];

    const auto edge_len = edge.length;

    float weight = fabs(c * (v2.p - v1.p)) / edge_len;

    weights[i] = weight;
  }
}

///////////////////////////////////////////////////////////////////////////////
void UnFlatTreeSplitter::assignWeightsImpl(model* m, vector<float>& weights,
    const Config& config) {
  FlatTreeSplitter::assignWeightsImpl(m, weights, config);

  for (auto& w : weights) {
    // edge needs to be avoided
    if (w <= 0)
      continue;
    w = 1.0 - w;
  }
}

//////////////////////////////////////////////////////////////////////
void SteepestEdgeSplitter::assignWeightsImpl(model* m, vector<float>& weights,
    const Config& config) {

  const auto c = this->genRandomUnitVector(config);

  int top_vertex_id = -1;

  {
    // find the top vertex w.r.t to random vector c
    float max_prod = -FLT_MAX;
    for (auto i = 0; i < m->v_size; ++i) {
      auto prod = c * Vector3d(m->vertices[i].p.get());
      if (prod > max_prod) {
        max_prod = prod;
        top_vertex_id = i;
      }
    }
  }

  unordered_set<int> cut_edges;

  for (auto i = 0; i < m->v_size; ++i) {

    // skip top vertex
    if (i == top_vertex_id)
      continue;

    const auto& v = m->vertices[i];

    // float vertex?
    if (v.m_e.empty())
      continue;

    float max_prod = -FLT_MAX;
    int steepest_eid = INT_MAX;
    int steepest_wid = INT_MAX;

    int steepest_eid2 = INT_MAX;
    int steepest_wid2 = INT_MAX;

    // find the steepest edge of this vertex
    for (auto eid : v.m_e) {
      const auto& e = m->edges[eid];

      // avoid diagonal edges
      if (config.less_cuts && (e.type == 'd' || e.type == 'p'))
        continue;

      const auto& wid = e.vid[0] == i ? e.vid[1] : e.vid[0];

      const auto& w = m->vertices[wid];

      const Vector3d vw = w.p - v.p;

      // prod in [-1, 1]
      float prod = c * vw / vw.norm();

      if (!m_steepest) {
        prod = 1.0 - fabs(prod);
      }

      if (prod > max_prod) {

        max_prod = prod;
        steepest_eid = eid;
        steepest_wid = wid;
      }

    }

    if (steepest_eid == INT_MAX) {
      cerr << "! Warning: steepest edge not found for vertex " << i << endl;
    } else {
      cut_edges.insert(steepest_eid);
    }
  }

// assign weight
  for (auto i = 0; i < m->e_size; ++i) {
    const auto& e = m->edges[i];

    if (cut_edges.count(i)) {
      weights[i] = EdgeWeight::CutEdge;
    } else {
      weights[i] = EdgeWeight::KeepEdge;
    }
  }
}

/////////////////////////////////////////////////////////////////////

void RandomSplitter::assignWeightsImpl(model* m, vector<float>& weights,
    const Config& config) {

  for (auto i = 0; i < m->e_size; i++) {
    const auto& edge = m->edges[i];
    const auto fid1 = edge.fid[0];
    const auto fid2 = edge.fid[1];

    // random weight with some offset
    float weight = mathtool::drand48() - EdgeWeight::DiagnalEdge;

    weights[i] = weight;

    cerr << weight << " ";
  }
}

/////////////////////////////////////////////////////////////////////

void RandomEdgeSplitter::assignWeightsImpl(model* m, vector<float>& weights,
    const Config& config) {

  for (auto i = 0; i < m->e_size; i++) {
    weights[i] = -EdgeWeight::DiagnalEdge; // add some offset
  }

  for (const auto& v : m->vertices) {
    // randomly select two edges for hyperbolic vertices
    if (v.hyperbolic) {
      vector<int> eids(v.m_e.begin(), v.m_e.end());
      int index0 = mathtool::drand48() * eids.size();
      int index1 = index0;
      while (index1 == index0) {
        index1 = mathtool::drand48() * eids.size();
      }
      int eid0 = eids[index0];
      int eid1 = eids[index1];
      weights[eid0] = EdgeWeight::CutEdge;
      weights[eid1] = EdgeWeight::CutEdge;
    }
  }

  for (const auto& v : m->vertices) {
    // randomly select one edge for normal vertex if no edge was selected
    if (!v.hyperbolic) {
      if (v.m_e.empty())
        continue;
      bool selected = false;
      for (uint eid : v.m_e)
        if (weights[eid] > 0) {
          selected = true;
          break;
        }

      if (selected)
        continue;

      vector<int> eids(v.m_e.begin(), v.m_e.end());
      int index0 = mathtool::drand48() * eids.size();
      int eid0 = eids[index0];
      weights[eid0] = EdgeWeight::CutEdge;
    }
  }
}

void BruteForceSplitter::init(int edges) {
  cout << "BruteForceSplitter::init edges = " << edges << endl;
  this->m_weights = vector<float>(edges);

  for (int i = 0; i < edges; ++i) {
    this->m_weights[i] = i;
  }
}

void BruteForceSplitter::assignWeightsImpl(model *m, vector<float>& weights,
    const Config& config) {
    if (!this->m_inited) {
      this->init(m->e_size);
      this->m_inited = true;
    } else {
      if (!std::next_permutation(this->m_weights.begin(),
          this->m_weights.end())) {
        cerr << "All possible permutation tried!" << endl;
        assert(false);
      }
    }
}

} /* namespace masc */
