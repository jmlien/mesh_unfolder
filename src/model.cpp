//------------------------------------------------------------------------------
//  Copyright 2007-2009 by Jyh-Ming Lien and George Mason University
//  See the file "LICENSE" for more information
//------------------------------------------------------------------------------

#include "model.h"

#include <string>
#include <algorithm>
using namespace std;

#include "ModelGraph.h"
#include "mathtool/Geometry.h"

#ifdef _WIN32
#pragma warning(disable:4244)
#endif

//-----------------------------------------------------------------------------
// two global scope models, P, Q
model P, Q;
model& getP() {
  return P;
}
model& getQ() {
  return Q;
}
//-----------------------------------------------------------------------------

//
// build the model from a file
// the file will contain an obj model
// once obj file is read, vertices and facets will be built
// other information associated with facets and vertices,
// as well as edges will be build also using model graph
//

bool model::build(const string & name, bool quiet)
{
  MeshReader* reader = nullptr;

  if (name.find(".obj") == name.length() - 4) {
    reader = new objReader();
  } else if (name.find(".off") == name.length() - 4) {
    reader = new offReader();
  }

  if (reader == nullptr) {
    cerr << "!Error! Unknown file type: " << name << endl;
    return false;
  }

  if (!reader->Read(name)) {
    delete reader;
    return false;
  }

  objModel& data = reader->getModel();
  this->texture_path = reader->GetTexturePath();
  this->material_path= reader->GetMaterialPath();

  bool result = build(data, quiet);

  delete reader;

  { //model name
    int start = name.find_last_of("/");
    if (start == string::npos)
      start = name.find_last_of("\\") + 1;
    if (start == string::npos)
      start = 0;
    int end = name.find_last_of(".");
    if (end == string::npos)
      end = name.length();
    this->name = name.substr(start, end - start);
  }

  return result;
}

bool model::build(const vector<Vector3d>& vertices,
                  const vector<vector<int>>& faces,
                  bool quiet)
{
  objModel data(vertices, faces);
  return this->build(data, quiet);
}

bool model::build(const objModel& data, bool quiet) {

  //allocate memory
  v_size = data.pts.size();
  t_size = data.polys.size();
  vt_size = data.textures.size();

  texture_pts.resize(vt_size);

  vertices.resize(v_size);   //
  tris.resize(t_size);     //
  //  assert(vertices && tris);        //make sure enough memory

  //copy vertices
  for (uint i = 0; i < v_size; i++) {
    vertices[i].p.set(&data.pts[i].x);
    vertices[i].bk_p = vertices[i].p;  //backup for modification
  }

  //cout << v_size << " vertices copied!" << endl;

  // copy texture vertices
  for (uint i = 0; i < vt_size; ++i) {
    texture_pts[i].set(data.textures[i].x, data.textures[i].y);
  }

  //cout << vt_size << " texture vertices copied!" << endl;

  //copy triangles
  int tid = 0;
  for (const polygon& poly : data.polys) {
    const auto& ids = poly.pts;
    //check if triangle
    if (ids.size() != 3) {
      cerr << "! Error: polygon " << tid << " is not a triangle, edge size is "
          << ids.size() << "." << endl;
      return false;
    }
    int vid = 0;
    for (auto j = ids.begin(); j != ids.end(); j++) {
      assert(*j >= 0 && *j < v_size);
      tris[tid].v[vid++] = *j;
      vertices[*j].m_f.push_back(tid);
    }

    // Has texture, copy texture vertex ids
    if (poly.textures.size() == 3) {
      for (int i = 0; i < 3; ++i) {
        tris[tid].vt[i] = poly.textures[i];
      }
    }

    tid++;
  }

  //cout << t_size << " triangles copied!" << endl;

  {  //build model graph and collect informations
    this->e_boundary_size = 0;

    CModelGraph G;
    G.doInit(this);
    //create edges
    e_size = G.getEdgeSize();
    CModelEdge * ptrE = G.getEdges();
    edges.resize(e_size);
//    assert(edges);
    for (uint i = 0; i < e_size; i++, ptrE = ptrE->getNext()) {
      int v1 = edges[i].vid[0] = ptrE->getStartPt();
      int v2 = edges[i].vid[1] = ptrE->getEndPt();

      const vector<int>&tmp_fid = ptrE->getFacets();

      edges[i].fid.insert(edges[i].fid.end(), tmp_fid.begin(), tmp_fid.end());

      if (tmp_fid.size() < 2) { //check if boundary edge
        //edges[i].fid[1]=edges[i].fid[0];
        edges[i].fid.push_back(edges[i].fid[0]);
        edges[i].type = 'b';		//bd
        this->e_boundary_size++;
      }

      //compute parallel vector
      edges[i].v = edges[i].bk_v = ptrE->getV();
      //edges[i].in_n[0]=ptrE->getInNormal(0); edges[i].bk_in_n[0]=edges[i].in_n[0];
      //edges[i].in_n[1]=ptrE->getInNormal(1); edges[i].bk_in_n[1]=edges[i].in_n[1];
      vertices[v1].m_e.push_back(i);
      vertices[v2].m_e.push_back(i);
    }            //end i

//facets
    vector<CModelFacet>& facets = G.getFacets();
    for (uint i = 0; i < t_size; i++) {
      tris[i].n = tris[i].bk_n = facets[i].n;
      for (int j = 0; j < 3; j++) {
        tris[i].e[j] = facets[i].m_Edge[j]->getID();
      }            //end j

      // compute face area
      double ss = 0.0;
      double s[3]; // edge length
      for (int j = 1; j <= 3; ++j) {
        s[j - 1] =
            (vertices[tris[i].v[j % 3]].p - vertices[tris[i].v[j - 1]].p).norm();
        ss += s[j - 1];
        tris[i].center += Vector3d(vertices[tris[i].v[j % 3]].p.get());
      } // end j

      // store center of the face
      tris[i].center = tris[i].center / 3.0;

      ss /= 2.0;

      tris[i].area = ss;
      for (int j = 0; j < 3; ++j) {
        tris[i].area *= (ss - s[j]);
      } // end j

      tris[i].area = std::sqrt(tris[i].area);

      this->surface_area += tris[i].area;
    }            //end i

//edge type
    ptrE = G.getEdges();
    for (uint i = 0; i < e_size; i++, ptrE = ptrE->getNext()) {

      // measure edge length
      const Point3d& p0 = this->vertices[edges[i].vid[0]].p;
      const Point3d& p1 = this->vertices[edges[i].vid[1]].p;

      const float length = (p1 - p0).norm();

      edges[i].length = length;

      edge& e = edges[i];

      this->total_edge_length += e.length;

      if (e.type == 'b')
        continue; //already know
      Vector3d& n1 = tris[edges[i].fid[0]].n;
      Vector3d& n2 = tris[edges[i].fid[1]].n;
      double d = n1 * n2;
      if (fabs(1 - d) < 1e-6) {
        e.type = 'p'; //plane
        e.folding_angle = 0;
      } else {
        double angle = acos(d);
        Vector3d vec = (n1 % n2).normalize();
        if (vec * e.v > 0) {
          e.type = 'c'; //convex
          e.folding_angle = angle;
        } else {
          e.type = 'r'; //reflex
          e.folding_angle = -angle;
        }
      }

      if (e.type != 'p')
        continue;



      // check whether the edge is diagonal edge of a quad
      const auto& t1 = tris[e.fid[0]];

      const auto& p11 = vertices[t1.v[0]].p;
      const auto& p12 = vertices[t1.v[1]].p;
      const auto& p13 = vertices[t1.v[2]].p;

      if (!mathtool::isRightTriangle(p11, p12, p13))
        continue;

      const auto& t2 = tris[e.fid[1]];

      const auto& p21 = vertices[t2.v[0]].p;
      const auto& p22 = vertices[t2.v[1]].p;
      const auto& p23 = vertices[t2.v[2]].p;

      if (!mathtool::isRightTriangle(p21, p22, p23))
        continue;

      const auto e_p1 = edges[t1.getPrevEdge(i)].v.normalize();
      const auto e_n1 = edges[t1.getNextEdge(i)].v.normalize();
      const auto e_p2 = edges[t2.getPrevEdge(i)].v.normalize();
      const auto e_n2 = edges[t2.getNextEdge(i)].v.normalize();

      double area1 = mathtool::triangleArea(p11, p12, p13);
      double area2 = mathtool::triangleArea(p21, p22, p23);

      const auto dot1 = e_p1 * e_n2;
      const auto dot2 = e_p2 * e_n1;



#define SMALL_NUMBER (1e-3)

      if (fabs(area1 - area2) / area1 < SMALL_NUMBER
          && fabs(dot1) < SMALL_NUMBER && fabs(dot2) < SMALL_NUMBER
          && e.type == 'p') {
        e.type = 'd';  // diagonal
        e.diagonal = true;
      }

    }

//vertex type
//    typedef list<uint>::iterator IT;
    for (uint i = 0; i < v_size; i++) {
      //int convex_e=-1;
      vertex& v = vertices[i];

      Vector3d evec;
      for (const auto eid : v.m_e) {
        edge& e = edges[eid];
        Vector3d vec = e.v;
        if (e.vid[1] == i)
          vec = -vec;
        evec = evec + vec;
      } //end ir

      Vector3d fvec;
      for (const auto fid : v.m_f) {
        triangle& t = tris[fid];
        fvec = fvec + t.n;
      } //end ir

      if (evec * fvec > 0)
        v.concave = true;

      double sum_angles = 0;

      for (const auto fid : v.m_f) {
        const auto& f = this->tris[fid];

        int vid1 = -1;
        int vid2 = -1;

        for (auto vid : f.v) {
          if (vid == i)
            continue;
          if (vid1 == -1)
            vid1 = vid;
          else
            vid2 = vid;
        }

        Vector3d e1 =
            (this->vertices[vid1].p - this->vertices[i].p).normalize();
        Vector3d e2 =
            (this->vertices[vid2].p - this->vertices[i].p).normalize();

        double angle = acos(e1 * e2);

        sum_angles += angle;
      }

      if (sum_angles > 2 * PI) {
        v.hyperbolic = true;
      }

    } //end i
  }

  if(!quiet)
    cout << "- v_size (vertex) = " << v_size << " t_size (triangle) = " << t_size << " e_size (edge) = "
         << e_size << " vt_size (uv) = "<<vt_size<<endl;

  this->compute_COM_R();

  if(!quiet) cout << "- COM = " << this->COM << " R = " << this->R << endl;

  return true;
}

//
//
//
// Transformation operations
//
// Rotation
// Negation
// Reverse (turn inside out)
//

void model::perturb(double noise) {

  //build the matrix
  Vector3d r(noise * mathtool::drand48(), noise * mathtool::drand48(),
      noise * mathtool::drand48());
  for (int i = 0; i < 3; i++) {
    if (mathtool::drand48() > 0.5)
      r[i] = -r[i];
  }

  Quaternion q(r.get());
  Matrix3x3 M = q.getMatrix();

  //rotate edges
  for (uint i = 0; i < e_size; i++) {
//edges[i].in_n[0]=(M*edges[i].in_n[0]).normalize();
//edges[i].in_n[1]=(M*edges[i].in_n[1]).normalize();
    edges[i].v = (M * edges[i].v).normalize();
  }
  //rotate facets
  for (uint i = 0; i < t_size; i++)
    tris[i].n = (M * tris[i].n).normalize();
}

void model::unperturb() {
  //build the matrx
  Matrix3x3 M(current_rot[0][0], current_rot[0][1], current_rot[0][2],
      current_rot[1][0], current_rot[1][1], current_rot[1][2],
      current_rot[2][0], current_rot[2][1], current_rot[2][2]);

  //rotate edges
  for (uint i = 0; i < e_size; i++) {
//edges[i].in_n[0]=(M*edges[i].bk_in_n[0]).normalize();
//edges[i].in_n[1]=(M*edges[i].bk_in_n[1]).normalize();
    edges[i].v = (M * edges[i].bk_v).normalize();
  }
  //rotate facets
  for (uint i = 0; i < t_size; i++)
    tris[i].n = (M * tris[i].bk_n).normalize();
}

void model::rotate(const Matrix2x2& M) {
  Vector2d tmp;

  //rotate vertices
  for (uint i = 0; i < v_size; i++) {
    tmp.set(vertices[i].bk_p[0], vertices[i].bk_p[1]);
    tmp = M * tmp;
    vertices[i].p.set(tmp[0], tmp[1]);
  }

  //rotate edges
  for (uint i = 0; i < e_size; i++) {
//for(int j=0;j<2;j++){
//    tmp.set(edges[i].bk_in_n[j][0],edges[i].bk_in_n[j][1]);
//    tmp=M*tmp;
//    edges[i].in_n[j].set(tmp[0],tmp[1]);
//}

    tmp.set(edges[i].v[0], edges[i].v[1]);
    tmp = M * tmp;
    edges[i].v.set(tmp[0], tmp[1]);
  }

  //rotate facets
  for (uint i = 0; i < t_size; i++) {
    tmp.set(tris[i].n[0], tris[i].n[1]);
    tmp = M * tmp;
    tris[i].n.set(tmp[0], tmp[1]);
  }
}

void model::rotate(const Matrix3x3& M) {
  Vector3d tmp;

  //rotate vertices
  for (uint i = 0; i < v_size; i++) {
    tmp.set(vertices[i].bk_p.get());
    vertices[i].p = M * tmp;
  }

  //rotate edges
  for (uint i = 0; i < e_size; i++) {
//edges[i].in_n[0]=(M*edges[i].bk_in_n[0]).normalize();
//edges[i].in_n[1]=(M*edges[i].bk_in_n[1]).normalize();
    edges[i].v = (M * edges[i].bk_v).normalize();
  }
  //rotate facets
  for (uint i = 0; i < t_size; i++)
    tris[i].n = (M * tris[i].bk_n).normalize();
}

void model::scale(double s) {
  //rotate vertices
  Point3d tmp;
  for (uint i = 0; i < v_size; i++) {
    tmp.set(vertices[i].bk_p.get());
    vertices[i].p.set(tmp[0] * s, tmp[1] * s, tmp[2] * s);
    vertices[i].bk_p.set(tmp[0] * s, tmp[1] * s, tmp[2] * s);
  }
}

void model::negate() {

  for (uint i = 0; i < v_size; i++) {
    Point3d& pt = vertices[i].p;
    pt.set(-pt[0], -pt[1], -pt[2]);

    Point3d& pt2 = vertices[i].bk_p;
    pt2.set(-pt2[0], -pt2[1], -pt2[2]);
  }

  for (uint i = 0; i < t_size; i++) {
    tris[i].n = -tris[i].n;
    tris[i].bk_n = -tris[i].bk_n;
    swap(tris[i].v[1], tris[i].v[2]);
//swap(tris[i].e[1],tris[i].e[2]);
  }

  for (uint i = 0; i < e_size; i++) {
    edge& e = edges[i];
    e.v = -e.v;
//e.in_n[0]=-e.in_n[0];
//e.in_n[1]=-e.in_n[1];
    e.bk_v = -e.bk_v;
//e.bk_in_n[0]=-e.bk_in_n[0];
//e.bk_in_n[1]=-e.bk_in_n[1];
  }
}

void model::reverse() {
  for (uint i = 0; i < t_size; i++) {
    tris[i].n = -tris[i].n;
    tris[i].bk_n = -tris[i].bk_n;
    swap(tris[i].v[1], tris[i].v[2]);
  }

  for (uint i = 0; i < e_size; i++) {
    edge& e = edges[i];
    if (e.type == 'r')
      e.type = 'c';
    else if (e.type == 'c')
      e.type = 'r';
  }
}

//get the neighbors of triangle t
/*
 void model::get_neighbors(triangle * t, list<triangle *>& nei)
 {
 for(int i=0;i<3;i++){
 edge& e=edges[t->e[i]];
 uint o_t=(&(tris[e.fid[0]])==f)?e.fid[1]:e.fid[0];
 if(&tris[o_t]==t) continue;
 nei.push_back(&tris[o_t]);
 }
 }
 */

void model::compute_COM_R() {
  this->COM = Vector3d(0, 0, 0);

  for (int i = 0; i < this->v_size; ++i) {
    this->COM = this->COM + Vector3d(this->vertices[i].p.get());
  }

  this->COM = this->COM / this->v_size;

  double r = 0;

  for (int i = 0; i < this->v_size; ++i) {
    double d = (Vector3d(this->vertices[i].p.get()) - this->COM).normsqr();

    if (d > r)
      r = d;
  }

  this->R = std::sqrt(r);

}

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

void model::getQuadIds(int eid, vector<int>* eids, vector<int>* vids) const {
  const edge& e = edges[eid];

  assert(e.type == 'd');

  const triangle& f1 = tris[e.fid[0]];
  const triangle& f2 = tris[e.fid[1]];

  // Collect edges of the triangles
  vector<uint> vs1(f1.v, f1.v + 3);
  vector<uint> vs2(f2.v, f2.v + 3);

  // Diagonal edge is v0 <-> v2

  // rotate to: v0 -> v1 -> v2
  while (vs1[1] == e.vid[0] || vs1[1] == e.vid[1]) {
    vs1.push_back(vs1.front());
    vs1.erase(vs1.begin());
  }

  // rotate to: v0 -> v2 -> v3
  while (vs2[2] == e.vid[0] || vs2[2] == e.vid[1]) {
    vs2.push_back(vs2.front());
    vs2.erase(vs2.begin());
  }

  // v0->v1->v2->v3
  vids->clear();
  for (const int vid : vs1)
    vids->push_back(vid);
  vids->push_back(vs2.back());

  const int e01 = getEdgeId(vids->at(0), vids->at(1));
  const int e12 = getEdgeId(vids->at(1), vids->at(2));
  const int e23 = getEdgeId(vids->at(2), vids->at(3));
  const int e30 = getEdgeId(vids->at(3), vids->at(0));

  *eids = vector<int> { e01, e12, e23, e30 };
}

/*

 v0    v3        v0    v3
 *-----*         *-----*
 | \ f2|         |f2 / |
 |  \  |   -->   |  /  |
 |f1 \ |         | /f1 |
 *-----*         *-----*
 v1    v2        v1    v2

 */
void model::flipEdge(int eid) {
  edge& e = edges[eid];

  // Only support flip diagonal edges.
  assert(e.type == 'd');

  vector<int> vids;
  vector<int> eids;

  this->getQuadIds(eid, &eids, &vids);

  int fid1 = e.fid[0];
  int fid2 = e.fid[1];

  triangle* f1 = &tris[fid1];
  triangle* f2 = &tris[fid2];

  // v1 should be in f1
  if (std::find(f1->v, f1->v + 3, vids[1]) == f1->v + 3) {
    std::swap(f1, f2);
    std::swap(fid1, fid2);
  }

  // Update Faces
  f1->v[0] = vids[1];
  f1->v[1] = vids[2];
  f1->v[2] = vids[3];

  f2->v[0] = vids[3];
  f2->v[1] = vids[0];
  f2->v[2] = vids[1];

  f1->e[0] = eids[1];
  f1->e[1] = eids[2];
  f1->e[2] = eid;

  f2->e[0] = eids[3];
  f2->e[1] = eids[0];
  f2->e[2] = eid;

  // Update vertices
  vertex* v0 = &vertices[vids[0]];
  vertex* v1 = &vertices[vids[1]];
  vertex* v2 = &vertices[vids[2]];
  vertex* v3 = &vertices[vids[3]];

  v0->m_e.remove(eid);
  v2->m_e.remove(eid);
  v1->m_e.push_back(eid);
  v3->m_e.push_back(eid);

  v0->m_f.remove(fid1);
  v2->m_f.remove(fid2);
  v1->m_f.push_back(fid2);
  v3->m_f.push_back(fid1);

  // Update edges
  edge* e01 = &edges[eids[0]];
  edge* e12 = &edges[eids[1]];
  edge* e23 = &edges[eids[2]];
  edge* e30 = &edges[eids[3]];

  e01->fid.erase(std::find(e01->fid.begin(), e01->fid.end(), fid1));
  e01->fid.push_back(fid2);

  e23->fid.erase(std::find(e23->fid.begin(), e23->fid.end(), fid2));
  e23->fid.push_back(fid1);

  // copy for modification
  edge new_edge = e;

  new_edge.vid[0] = vids[1];
  new_edge.vid[1] = vids[3];
  edges[eid] = new_edge;

}

void model::flipEdge(int vid1, int vid2) {
  int eid = getEdgeId(vid1, vid2);
  assert(eid >= 0);
  return flipEdge(eid);
}

int model::getEdgeId(int vid1, int vid2) const {
  for (const auto eid : this->vertices[vid1].m_e) {
    const edge& e = this->edges[eid];
    if (e.vid[0] == vid2 || e.vid[1] == vid2)
      return eid;
  }
  return -1;
}

int model::getFaceId(int eid1, int eid2) const {
  const edge& e1 = this->edges[eid1];

  for (int fid : e1.fid) {
    const triangle& tri = this->tris[fid];
    for (int eid : tri.e) {
      if (eid == eid2)
        return fid;
    }
  }

  return -1;
}

bool model::isBorderVertex(int vid) const {
  const vertex& v = this->vertices[vid];
  for (const int eid : v.m_e) {
    const edge& e = this->edges[eid];
    if (e.type == 'b')
      return true;
  }
  return false;
}

bool model::isNeighbor(int eid1, int eid2) const {
  const edge& e1 = this->edges[eid1];
  const edge& e2 = this->edges[eid2];

  for (int fid : e1.fid) {
    const triangle& tri = this->tris[fid];
    for (int eid : tri.e) {
      if (eid == eid2)
        return true;
    }
  }

  for (int fid : e2.fid) {
    const triangle& tri = this->tris[fid];
    for (int eid : tri.e) {
      if (eid == eid1)
        return true;
    }
  }

  return false;
}

bool model::isNeighborFaces(int fid1, int fid2) const {
  const triangle& f1 = this->tris[fid1];
  const triangle& f2 = this->tris[fid2];

  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      if (f1.e[i] == f2.e[j])
        return true;
    }
  }

  return false;
}

bool model::getHalfRing(int eid, int fid, int vid, vector<int>* eids) const {

  const edge& start_edge = this->edges[eid];
  const vertex& vertex = this->vertices[vid];
  const triangle& face = this->tris[fid];

// start edge should not be a boundary edge.
  assert(start_edge.type != 'b');
// vid must be on the boundary
  assert(isBorderVertex(vid));

  eids->clear();
  eids->push_back(eid);

  set<int> visted = { eid };

  int last_eid = eid;
  bool found = true;

  while (found) {
    found = false;
    for (int cur_eid : vertex.m_e) {

      if (!visted.count(cur_eid) && isNeighbor(cur_eid, last_eid)) {

        // the second edge in half ring should be in face
        if (visted.size() == 1) {
          const edge& cur_edge = this->edges[cur_eid];
          if (cur_edge.fid[0] != fid && cur_edge.fid[1] != fid)
            continue;
        }

        // Move to the next edge.
        last_eid = cur_eid;
        eids->push_back(last_eid);
        visted.insert(last_eid);
        found = true;
        break;

      }
    }
  }

  const edge& end_edge = this->edges[last_eid];
  // last edge should be on the boundary
  assert(end_edge.type == 'b');

  return true;
}

int model::getEdgeIdByFids(int fid1, int fid2) const
{
  const triangle& f1 = this->tris[fid1];
  const triangle& f2 = this->tris[fid2];
  for (int eid : f1.e) {
    const edge& e = this->edges[eid];
    if (std::find(e.fid.begin(), e.fid.end(), (unsigned int) fid2)
        != e.fid.end()) {
      return eid;
    }
  }

  for (int eid : f2.e) {
    const edge& e = this->edges[eid];
    if (std::find(e.fid.begin(), e.fid.end(), (unsigned int) fid1)
        != e.fid.end()) {
      return eid;
    }
  }
  return -1;
}

void model::updateHalfRing(const vector<int>& ring_eids, int old_vid,
    int new_vid) {

  vertex* old_vertex = &this->vertices[old_vid];
  vertex* new_vertex = &this->vertices[new_vid];

  // Use the new vertex on the half ring
  int last_r_eid = ring_eids.front();
  for (int i = 1; i < ring_eids.size(); ++i) {
    int r_eid = ring_eids[i];
    edge& r_e = this->edges[r_eid];

    // Set the new vid for each edge in the half ring
    if (r_e.vid[0] == new_vid || r_e.vid[1] == new_vid) {
      // should not happen!!!
      assert(false);
    } else if (r_e.vid[0] == old_vid) {
      r_e.vid[0] = new_vid;
    } else {
      r_e.vid[1] = new_vid;
    }

    //get the triangle incident to both last_r_eid, r_eid
    int r_fid = this->getFaceId(last_r_eid, r_eid);

    assert(r_fid>=0);


    triangle& r_f = this->tris[r_fid];

    bool replaced = false;

    // Set the new vid for each face in the half ring
    for (int j = 0; j < 3; ++j) {
      if (r_f.v[j] == old_vid) {
        r_f.v[j] = new_vid;
        replaced = true;
        break;
      }
    }

    assert(replaced);

    old_vertex->m_f.remove(r_fid);
    new_vertex->m_f.push_back(r_fid);


    old_vertex->m_e.remove(r_eid);
    new_vertex->m_e.push_back(r_eid);

    last_r_eid = r_eid;
  }
}

void model::cutEdge(int eid_to_cut) {

  const edge& edge_to_cut = this->edges[eid_to_cut];
  int vid1 = edge_to_cut.vid[0];
  int vid2 = edge_to_cut.vid[1];

  assert(vid1 != vid2);

  if (!this->isBorderVertex(vid1) && !this->isBorderVertex(vid2)) {
    // case 1, both vertices are non-border vertices
    // add an new edge in opposite direction.

    int old_eid = this->getEdgeId(vid1, vid2);
    edge& old_edge = this->edges[old_eid];

    // Should have two faces...
    assert(old_edge.fid.size() == 2);

    uint fid1 = old_edge.fid[0];
    uint fid2 = old_edge.fid[1];

    assert(fid1 != fid2);

    triangle& f1 = this->tris[fid1];

    // f2 will use the new edge
    triangle& f2 = this->tris[fid2];

    // set folding angle to 0 for cut edges.
    old_edge.folding_angle = 0.0;

    // copy every thing form the old edge
    edge new_edge = old_edge;
    new_edge.type = 'b'; // new edge is a border edge
    new_edge.fid = {fid2};
    new_edge.parent_id = old_eid;

    old_edge.type = 'b'; // old edge becomes a border edge
    old_edge.fid = {fid1};

    const uint new_edge_id = this->edges.size();

    // update f2's edge list
    for (int i = 0; i < 3; ++i)
      if (f2.e[i] == old_eid)
        f2.e[i] = new_edge_id;

    // update v1,v2's edge list
    vertex& v1 = vertices[vid1];
    vertex& v2 = vertices[vid2];

    v1.m_e.push_back(new_edge_id);
    v2.m_e.push_back(new_edge_id);

    this->edges.push_back(new_edge);
    ++this->e_size;

  } else if (this->isBorderVertex(vid1) ^ this->isBorderVertex(vid2)) {
    // case 2, one vertex is on the boundary
    // one vertex is on the boundary

    // Get the edge to cut
    int old_eid = eid_to_cut;
    edge& old_edge = this->edges[old_eid];

    // Get the boundary vertex and the inner vertex
    int b_vid = this->isBorderVertex(vid1) ? vid1 : vid2;
    int i_vid = (b_vid == vid1) ? vid2 : vid1;

    // boundary vertex
    vertex* b_vertex = &this->vertices[b_vid];

    // inner vertex
    vertex* i_vertex = &this->vertices[i_vid];

    // Get the half ring
    // All faces, will use the new vertex
    vector<int> eids;
    this->getHalfRing(old_eid, (int) old_edge.fid[1], b_vid, &eids);
    // At least two edges
    assert(eids.size() >= 2);

    int fid2 = old_edge.fid[1]; // the face on the half ring.
    int fid1 = old_edge.fid[0]; // the opposite face on old_edge.

    // Create a new vertex
    int new_vid = this->v_size;
    vertex new_vertex;
    new_vertex.p = b_vertex->p; // Same position.
    new_vertex.parent_id = b_vid;

    // Create a new edge
    int new_eid = this->e_size;
    edge new_edge = old_edge;
    new_edge.type = 'b';
    new_edge.vid[0] = (old_edge.vid[0] == i_vid ? i_vid : new_vid);
    new_edge.vid[1] = (old_edge.vid[1] == i_vid ? i_vid : new_vid);
    new_edge.fid = {(uint)fid2};
    new_edge.parent_id = old_eid;

    old_edge.type = 'b'; // set old edge as boundary edge.
    old_edge.fid = {(uint)fid1}; // remove one face.

    new_vertex.m_e.push_back(new_eid);
    i_vertex->m_e.push_back(new_eid);

    // Add to the model.
    this->edges.push_back(new_edge);
    this->vertices.push_back(new_vertex);

    // reassign after update the push_back operation...
    b_vertex = &this->vertices[b_vid];
    i_vertex = &this->vertices[i_vid];

    // new vertex
    vertex* n_vertex = &(this->vertices.back());

    ++this->v_size;
    ++this->e_size;

    // Set f2 with the new edge
    triangle& f2 = this->tris[fid2];
    for (int i = 0; i < 3; ++i) {
      if (f2.e[i] == old_eid)
        f2.e[i] = new_eid;
    }

    // Use the new vertex on the half ring
    int last_r_eid = new_eid;
    for (int i = 1; i < eids.size(); ++i) {
      int r_eid = eids[i];
      edge& r_e = this->edges[r_eid];

      assert(r_e.vid[0] == b_vid || r_e.vid[1] == b_vid);

      // Set the new vid for each edge in the half ring
      if (r_e.vid[0] == new_vid || r_e.vid[1] == new_vid) {
        // should not happen!!!
        assert(false);
      } else if (r_e.vid[0] == b_vid) {
        r_e.vid[0] = new_vid;
      } else {
        r_e.vid[1] = new_vid;
      }

      int r_fid = this->getFaceId(last_r_eid, r_eid);

      assert(r_fid >= 0);


      triangle& r_f = this->tris[r_fid];

      // Set the new vid for each face in the half ring
      for (int j = 0; j < 3; ++j) {
        if (r_f.v[j] == b_vid)
          r_f.v[j] = new_vid;
      }

      b_vertex->m_f.remove(r_fid);
      n_vertex->m_f.push_back(r_fid);


      b_vertex->m_e.remove(r_eid);
      n_vertex->m_e.push_back(r_eid);

      last_r_eid = r_eid;
    }

  } else {
    // case 3, both vertices are on the boundary, this should only happens when genus > 0 and for exact 2*genus times.
    //         or when the model is open

    const int old_eid = eid_to_cut;
    edge* old_edge = &this->edges[old_eid];

    assert(old_edge->type != 'b');

    int fid1 = old_edge->fid[0];
    int fid2 = old_edge->fid[1];

    int vid1 = old_edge->vid[0];
    int vid2 = old_edge->vid[1];

    vector<int> ring_eids1;
    vector<int> ring_eids2;

    // Get half rings on the f2' side
    this->getHalfRing(old_eid, fid2, vid1, &ring_eids1);
    this->getHalfRing(old_eid, fid2, vid2, &ring_eids2);

    int vid1p = this->v_size;
    int vid2p = this->v_size + 1;

    // create new vertices
    vertex v1p;
    vertex v2p;

    v1p.p = this->vertices[vid1].p;
    v2p.p = this->vertices[vid2].p;
    v1p.parent_id = vid1;
    v2p.parent_id = vid2;

    // create a new edge
    int new_eid = this->e_size;
    edge new_edge = *old_edge; // copy from old edge
    new_edge.type = 'b'; // become a boundary edge
    new_edge.vid[0] = vid1p;
    new_edge.vid[1] = vid2p;
    new_edge.fid = {(uint)fid2};
    new_edge.parent_id = old_eid;

    // Add new vertices to the model.
    this->vertices.push_back(v1p);
    this->vertices.push_back(v2p);

    // Add new edge to the model.
    this->edges.push_back(new_edge);

    this->v_size += 2;
    this->e_size += 1;

    // Init new vertices
    vertex* ptr_v1p = &this->vertices[vid1p];
    vertex* ptr_v2p = &this->vertices[vid2p];

    ptr_v1p->m_e.push_back(new_eid);
    ptr_v2p->m_e.push_back(new_eid);
    ptr_v1p->m_f.push_back(fid2);
    ptr_v2p->m_f.push_back(fid2);

    // Use the new vertex on the half ring
    this->updateHalfRing(ring_eids1, vid1, vid1p);
    this->updateHalfRing(ring_eids2, vid2, vid2p);

    // Update f2 using the new edge. (new vertices should be updated in updateHalfRing)
    triangle& f2 = this->tris[fid2];
    for(int i=0;i<3;++i) {
      if(f2.e[i] == old_eid) f2.e[i] = new_eid;
    }

    // update old edge, should be the last step
    old_edge = &this->edges[old_eid];
    old_edge->type = 'b'; // also become a boundary edge
    old_edge->fid = {(uint)fid1};
  }

}

bool model::isEdgeCCW(int parent_fid, int vid1, int vid2) const {
  for (auto k = 0; k < 3; k++) {
    const auto fvid1 = this->tris[parent_fid].v[k];
    const auto fvid2 = this->tris[parent_fid].v[(k + 1) % 3];

    if (fvid1 == vid1 && fvid2 == vid2)
      return true;
  }

  return false;
}

void model::buildDualGraph() {

  // Return if already built.
  if (!dual_graph_weights.empty())
    return;

  const int INF_WEIGHTS = 9999;

  const int SQUARES = this->t_size / 2;

  this->dual_graph_weights = vector<vector<int>>(SQUARES,
      vector<int>(SQUARES, INF_WEIGHTS));

  for (int i = 0; i < dual_graph_weights.size(); ++i) {
    dual_graph_weights[i][i] = 0;
  }

  // dual node id
  int nid = 0;

  this->nid2eid = vector<int>(SQUARES, -1);
  this->eid2nid = vector<int>(e_size, -1);

  // face id -> edge id
  vector<int> f2e(this->t_size, -1);

  // First pass, collect mapping.
  for (int eid = 0; eid < e_size; ++eid) {
    const edge& e = edges[eid];
    if (e.type != 'd')
      continue;
    f2e[e.fid[0]] = eid;
    f2e[e.fid[1]] = eid;
    eid2nid[eid] = nid;
    nid2eid[nid] = eid;

    ++nid;
  }

  // Second pass, assign weights.
  for (int eid = 0; eid < e_size; ++eid) {
    const edge& e = edges[eid];
    if (e.type == 'b' || e.type == 'd')
      continue;
    int nid1 = eid2nid[f2e[e.fid[0]]];
    int nid2 = eid2nid[f2e[e.fid[1]]];

    // using uniform weight
    int weight = 1;

    dual_graph_weights[nid1][nid2] = weight;
    dual_graph_weights[nid2][nid1] = weight;
  }
}

void model::printObj(ostream& out) const
{
  const double tiny = 1e-11;
  out<<"# created by Mesh Unfolder from masc.cs.gmu.edu\n"
     <<"# there are "<<vertices.size()<<" vertices and "<<tris.size()<<" triangles\n"
     <<"mtllib "<<material_path<<"\n";

  for (const auto& v : vertices) {
    out << "v";
    out << " " << ((fabs(v.p[0]) < tiny) ? 0 : v.p[0]);
    out << " " << ((fabs(v.p[1]) < tiny) ? 0 : v.p[1]);
    out << " " << ((fabs(v.p[2]) < tiny) ? 0 : v.p[2]);
    out << "\n";
  }

  if(!texture_pts.empty())
  {
    for (const auto& v : texture_pts) {
      out << "vt "<<v<<"\n";
    }
  }

  if(!texture_pts.empty())
    for(const auto& f : tris) out << "f " << f.v[0] + 1 << "/" <<f.vt[0] + 1 << " " << f.v[1] + 1<<"/"<<f.vt[1] + 1 << " " << f.v[2] + 1<<"/"<<f.vt[2] + 1 << "\n";
  else for(const auto& f : tris) out << "f " << f.v[0] + 1 << " " << f.v[1] + 1 << " " << f.v[2] + 1 << "\n";
}

void model::saveObj(const string& path) const {
  ofstream out;
  out.open(path, ofstream::out);
  if (!out.good()) {
    cerr << "!Error! failed to open " << path << endl;
    return;
  }
  this->printObj(out);
  out.close();
  cerr << " - model output to " << path << endl;
}

model * model::create_submodel(const list<uint> & fids)
{
    objModel data;

    // UV coordinates
    bool has_uv= this->texture_pts.empty()==false;

    //collect all vertices
    unordered_map<uint, uint> vids;
    unordered_map<uint, uint> uvs;
    vids.reserve(this->v_size);
    uvs.reserve(this->v_size);

    for (auto fid : fids)
    {
        for (short d = 0; d<3; d++){
          vids[this->tris[fid].v[d]] = 0;
          if(has_uv) uvs[this->tris[fid].vt[d]] = 0;
        }
    }

    //get a list of points
    uint new_vid = 0;
    for (auto & vid : vids)
    {
        Vpt pt;
        auto & pos = this->vertices[vid.first].p;
        pt.x = pos[0];
        pt.y = pos[1];
        pt.z = pos[2];
        data.pts.push_back(pt);
        vid.second = new_vid++;
    }

    //go through the uvs
    uint new_uv=0;
    for (auto & id : uvs)
    {
      V uv;
      auto & pos=this->texture_pts[id.first];
      uv.x = pos[0];
      uv.y = pos[1];
      data.textures.push_back(uv);
      id.second = new_uv++;
    }

    //get a list of faces
    for (auto& fid : fids)
    {
        polygon poly;
        for (short d = 0; d < 3; d++){
            auto vid=vids[this->tris[fid].v[d]];
            poly.pts.push_back(vid);
            if(has_uv){
              auto uv=uvs[this->tris[fid].vt[d]];
              poly.textures.push_back(uv);
            }
        }
        data.polys.push_back(poly);
    }

    data.compute_v_normal();

    //build mesh
    model * subm = new model();
    if (subm->build(data, true) == false)
    {
        cerr << "! Error: Failed to build a model from face list" << endl;
        return NULL;
    }

    // Path of the texture file.
    subm->texture_path =this->texture_path;
    subm->material_path=this->material_path;

    //register
    auto fid_it = fids.begin();
    for (int fid = 0; fid < subm->t_size; fid++)
    {
        int sfid=*fid_it;
        int ssfid=tris[sfid].source_fid;
        subm->tris[fid].source_fid = (ssfid==-1)?sfid:ssfid;
        ++fid_it;
    }

    return subm;
}
