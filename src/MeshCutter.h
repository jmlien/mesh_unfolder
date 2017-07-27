/*
 * MeshCutter.h
 *
 *  Created on: Oct 14, 2016
 *      Author: zxi
 */

#ifndef SRC_MESHCUTTER_H_
#define SRC_MESHCUTTER_H_

#include <vector>
#include <utility>
#include <string>
using namespace std;

// forward declaration
struct model;

namespace masc {
namespace unfolding {
namespace mesh_cutter {
// Cut the mesh along given edges
// parameters:
//  mesh, the mesh to cut,
//  dual_edges: {<fid1, fid2>}
bool CutMesh(model* mesh, const vector<pair<int, int>>& dual_edges);

// Save vertex mapping list
// org_vid -> {org_vid, new_vid1, new_vid2, ... }
void SaveVertexMapping(const string& path, const model* model);
} // namespace mesh_cutter
} // namespace unfolding
} // namespace masc

#endif /* SRC_MESHCUTTER_H_ */
