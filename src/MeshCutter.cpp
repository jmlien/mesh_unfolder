/*
 * MeshCutter.cpp
 *
 *  Created on: Oct 14, 2016
 *      Author: zxi
 */

#include "MeshCutter.h"

#include <fstream>
#include <unordered_map>

#include "model.h" // for model

namespace masc {
namespace unfolding {
namespace mesh_cutter {

bool CutMesh(model* mesh, const vector<pair<int, int>>& dual_edges) {

  set<pair<int, int>> edges_to_cut(dual_edges.begin(), dual_edges.end());

  while (!edges_to_cut.empty()) {
    bool cutted = false;

    for (const auto& dual_edge : edges_to_cut) {
      int eid = mesh->getEdgeIdByFids(dual_edge.first, dual_edge.second);
      assert(eid >= 0);
      const edge& e = mesh->edges[eid];
      int vid1 = e.vid[0];
      int vid2 = e.vid[1];

      if (!mesh->isBorderVertex(vid1) && !mesh->isBorderVertex(vid2)) {
        // postpone the cut...
        continue;
      }

      mesh->cutEdge(eid);

      edges_to_cut.erase(dual_edge);

      cutted = true;

      break;
    }

    // Mesh is closed ....
    if (!cutted) {
      for (const auto& dual_edge : edges_to_cut) {
        int eid = mesh->getEdgeIdByFids(dual_edge.first, dual_edge.second);
        assert(eid >= 0);
        const edge& e = mesh->edges[eid];
        int vid1 = e.vid[0];
        int vid2 = e.vid[1];

        if (mesh->isBorderVertex(vid1) || mesh->isBorderVertex(vid2)) {
          // postpone the cut...
          continue;
        }

        mesh->cutEdge(eid);

        edges_to_cut.erase(dual_edge);

        cutted = true;

        break;
      }
    }

  }

  return true;
}

void SaveVertexMapping(const string& path, const model* model) {

  cerr << " - Saving vertex mapping to " << path << endl;

  ofstream out(path);

  unordered_map<int, vector<int>> vertex_mapping;

  for (size_t vid = 0; vid < model->v_size; ++vid) {
    const vertex& v = model->vertices[vid];
    int pid = vid;
    while (model->vertices[pid].parent_id != UINT_MAX) {
      pid = model->vertices[pid].parent_id;
    }
    vertex_mapping[pid].push_back(vid);
  }

  for (const auto& kv : vertex_mapping) {
    if (kv.second.size() < 2)
      continue;
    out << kv.first;
    for (const auto& vid : kv.second) {
      if (vid == kv.first)
        continue;
      out << " " << vid;
    }
    out << endl;
  }

  out.close();
}

} // namespace mesh_cutter
} // namespace unfolding
} // namespace masc
