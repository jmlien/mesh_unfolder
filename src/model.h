//------------------------------------------------------------------------------
//  Copyright 2007-2009 by Jyh-Ming Lien and George Mason University
//  See the file "LICENSE" for more information
//------------------------------------------------------------------------------

#ifndef _BF_MODEL_H_
#define _BF_MODEL_H_

#include <climits>
#include <cassert>

#include <string>
#include <set>
#include <unordered_map>
#include <utility>
#include <list>
using namespace std;

#include "mathtool/Point.h"
#include "mathtool/Vector.h"
#include "mathtool/Matrix.h"
#include "mathtool/Quaternion.h"
using namespace mathtool;

//#include "gmap.h"
#include "objReader.h"

typedef unsigned int uint;

//a triangle of the model
struct triangle {

  triangle() {
    cluster_id = 0;
    overlapped = false;
    parent_id = -1;
    source_fid = -1;
    path_len = 0.0;
    area = 0.0;
    score = 0.0;
  }

  uint v[3]; //vertex id
  uint e[3]; //edge id
  uint vt[3]; //texture vertex id;

  Vector3d n; //normal

  // source face id
  int source_fid;

  uint cluster_id; //id of the corresponding facet in kway-union.h

  //backups
  Vector3d bk_n;

  // for unfolding
  bool overlapped;

  // parent id, -1 means root face
  int parent_id;

  // path length from root
  double path_len;

  // area of the face
  double area;

  // center of mass
  Vector3d center;

  float score;

  int getPrevEdge(int eid) const {
    for (int i = 0; i < 3; ++i) {
      if (e[i] == eid)
        return e[(i + 2) % 3];
    }
    return -1;
  }

  int getNextEdge(int eid) const {
    for (int i = 0; i < 3; ++i) {
      if (e[i] == eid)
        return e[(i + 1) % 3];
    }
    return -1;
  }

};

//a vertex of the model
struct vertex {
  vertex() {
    concave = false;
    hyperbolic = false;
    score = 0.0;
    cut_src_id=UINT_MAX;
  }
  Point3d p;  //position
  list<uint> m_f;
  list<uint> m_e; //a list of edges

  //backups
  Point3d bk_p;

  //if concave, set to true
  bool concave;

  //whether the vertex is hyperbolic
  bool hyperbolic;

  // weighted score
  float score;

  //if this is a vertex created due to cut,
  //store where this vertex is from
  uint cut_src_id;

  //when this model is a sub-model created from a model
  //this source_vid points to the id of the original vertex in the original model
  uint source_vid;
};

//an edge of the model
struct edge {
  edge() {
    type = 'x';
    vid[0] = vid[1] = UINT_MAX;
    length = 0.0f;
    folding_angle = 0.0;
    cutted = true;
    diagonal = false;
    cut_twin_id=UINT_MAX;
    source_eid=UINT_MAX;
  }

  //given an incident face id (ofid)
  //return the other incident face id
  //only one id is returned if this is a non-manifold edge
  uint otherf(uint ofid) const
  {
    for(uint id : fid)
    {
      if(id!=ofid) return id;
    }
    return -1; //border edge
  }

  uint vid[2];
  vector<uint> fid; //incident face ids

  Vector3d v;       //parallel vector
  //Vector3d in_n[2]; //inface normals

  //backups
  Vector3d bk_v;       //parallel vector
  //Vector3d bk_in_n[2]; //inface normals

  //type, c-convex, r-reflex, b-border, p-plane, d-diagonal
  char type;

  // added for origami
  double folding_angle;
  float length;
  bool cutted;      // whether the edge is cut or not
  bool diagonal;    // whether the edge is a diagonal edge

  uint cut_twin_id;   //if this is a cut border edge, store where the twin of this cut edge

  uint source_eid; //if this edge is created from another model, this points to the original eid
};

struct model {

  //initialization
  model()
  {
    v_size = e_size = t_size = vt_size = e_boundary_size = 0;
    //spheres=NULL;
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
        current_rot[i][j] = 0;
    current_rot[0][0] = current_rot[1][1] = current_rot[2][2] = 1;
    R = 1.0;

    surface_area = 0.0;
    total_edge_length = 0.0;
    source=NULL;
  }

  ~model() {
  }

  void destroy() {
    tris.clear();
    edges.clear();
    vertices.clear();
    v_size = e_size = t_size = 0;
  }

  //build from model file: OBJ/OFF
  bool build(const std::string& name, bool quiet=false);

  //build from obj_model
  bool build(const masc::obj::objModel& obj_model, bool quiet=false);

  //build from vertices and faces
  bool build(const vector<Vector3d>& vertices,
  const vector<vector<int>>& faces, bool quiet=false);

  //create a submodel from a list of face ids
  model * create_submodel(const list<uint> & fids);

  //rotate points
  void rotate(const Matrix2x2& m);
  void rotate(const Matrix3x3& M);

  //scale the model
  void scale(double s);

  //negate point/facets ...
  void negate();

  //reverse facets ...
  void reverse();

  //perturb
  void perturb(double amount);
  void unperturb();

  // compute COM and R
  void compute_COM_R();

  // Get the eids/vids of a quad for a given the eid of diagonal edge.
  /*

    v0    v3        v3    v2
    *-----*         *-----*
    | \ f2|         |f2 / |
    |  \  |         |  /  |
    |f1 \ |         | /f1 |
    *-----*         *-----*
    v1    v2        v0    v1

    eids = {v0v1, v1v2, v2v3, v3v0}

  */
  void getQuadIds(int eid, vector<int>* eids, vector<int>* vids) const;

  // flip the edge, the edge to flip must be a diagonal edge.
  // modify the mesh in place, eid will not be changed.
  void flipEdge(int eid);
  void flipEdge(int vid1, int vid2);

  //void get_neighbors(triangle * t, list<triangle *>& nei);

  // get the edge id (0-based) for a given pair of vertices.
  // return -1 if no edge found.
  int getEdgeId(int vid1, int vid2) const;

  // get the edge id (0-based) for a given pair of faces.
  // return -1 if no edge found.
  int getEdgeIdByFids(int fid1, int fid2) const;

  // Get the face id (0-based) for a given pair of edges.
  int getFaceId(int eid1, int eid2) const;

  // cut an edge.
  void cutEdge(int eid);

  // Check whether a given edge (vid1->vid2) is counter clock wise with respsect to a given parent face (parent_fid)
  bool isEdgeCCW(int parent_fid, int vid1, int vid2) const;

  // Check whether two edges are neighbors or not.
  bool isNeighbor(int eid1, int eid2) const;

  // Check whether two faces are neighbors or not. e.g sharing an edge.
  bool isNeighborFaces(int fid1, int fid2) const;

  // Collect the edge ids in order around vid from the eid to any boundary edge by crossing fid
  // eids are ordered from eid to a boundary edge
  // vid must be on the boundary
  bool getHalfRing(int eid, int fid, int vid, vector<int>* eids) const;

  // Update the half ring with the new vertex.
  void updateHalfRing(const vector<int>& ring_eids, int old_vid, int new_vid);

  // Check whether a vertex is on the border of the mesh or not.
  bool isBorderVertex(int vid) const;

  // Print the model to out in obj format
  void printObj(ostream& out) const;

  // Save the model to an obj file.
  void saveObj(const string& path) const;

  //make this model manifold if it is non-manifold
  void makeManifold();

  //split a non-manifold vertex into multiple vertices
  void splitNonManifoldVertex(int old_vid);

  //data
  vector<vertex> vertices;  //vertices
  vector<triangle> tris;      //triangles
  vector<edge> edges;     //edges
  uint v_size;
  uint vt_size;     //texture coordinates size
  uint e_size;
  uint t_size;
  uint e_boundary_size; // size of boundary edges

  //current orientation
  double current_rot[3][3];

  Vector3d COM;
  double R;
  string name;

  // total surface area
  double surface_area;

  //
  // textures
  //

  // UV coordinates
  std::vector<Vector2d> texture_pts;

  // Path of the texture file.
  std::string texture_path;

  // Path to the material file
  std::string material_path;

  // total edge length
  double total_edge_length;

  /////////////////////////////////////////////////////////////
  // Should only used for dual graphs on squares meshes
  /////////////////////////////////////////////////////////////

  // build a dual graph
  void buildDualGraph();

  // fid1, fid2 -> weight, weights of dual graph
  vector<vector<int>> dual_graph_weights;

  // Diagonal edge id -> dual graph node id
  vector<int> eid2nid;

  // Dual graph node id -> Diagonal edge id
  vector<int> nid2eid;

  ////////////////////////////////////////////////////////////
  model * source; //if this is a sub-model of another model, source points where this model is from
};

// {<fid1,fid2>}
typedef vector<pair<int, int>> HamiltonianPath;
typedef unordered_map<uint, unordered_map<uint, float> > GRAPH; // adjacency list of weight <v_i, <v_j, w_ij>>
typedef vector<vector<pair<uint, Vector3d> > > MESH; // coordinates faces[f_i][k] = <vid, coordinates> (k = 0..2)

#endif //_BF_MODEL_H_
