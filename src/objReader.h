//------------------------------------------------------------------------------
//  Copyright 2007-2019 by Jyh-Ming Lien and George Mason University
//  See the file "LICENSE" for more information
//------------------------------------------------------------------------------

#ifndef _OBJ_READER_H_
#define _OBJ_READER_H_

#include <cctype>
#include <cmath>
#include <cassert>

#include <string>
#include <fstream>
#include <iostream>
#include <list>
#include <vector>

using namespace std;

#include "mathtool/Point.h"
#include "mathtool/Quaternion.h"
#include "mathtool/Vector.h"
using namespace mathtool;

#include "util/MTLReader.h"
namespace masc{
  namespace obj{

    class Vpt {
    public:
      Vpt() {
        x = y = z = nx = ny = nz = 0;
      }

      void normalize() {
        double norm = (double) sqrt(nx * nx + ny * ny + nz * nz);
        nx /= norm;
        ny /= norm;
        nz /= norm;
      }

      double x, y, z, nx, ny, nz;
    };

    class V {
    public:
      V() {
        x = y = z = 0;
      }
      double x, y, z;
    };

    istream& operator>>(istream& in, V& v);
    ostream& operator<<(ostream& out, const V& v);

    // Represent a face in the mesh
    class polygon {
    public:
      // Id of vertices
      vector<int> pts;

      // Id of vertex normals
      vector<int> normals;

      // Id of texture coordinates
      vector<int> textures;
    };

    class objModel {
    public:
      objModel() {
      }

      // build from vertices and faces
      objModel(const vector<Vector3d>& vertices, const vector<vector<int>>& faces) {
        Vpt pt;
        //  vertices
        for (const auto& v : vertices) {
          pt.x = v[0];
          pt.y = v[1];
          pt.z = v[2];
          this->pts.push_back(pt);
        }

        // faces
        for (const auto& f : faces) {
          polygon ply;
          ply.pts.push_back(f[0]);
          ply.pts.push_back(f[1]);
          ply.pts.push_back(f[2]);
          this->polys.push_back(ply);
        }

        this->compute_v_normal();
      }

      void compute_v_normal() {
        //check if normal information is valid from the obj file
        if (normals.empty()) { //compute normal
          for (polygon& i : polys) {
            //get 3 points, compute normal and assign to all vertices
            auto pi = i.pts.begin();
            vector<Point3d> v3;
            for (; pi != i.pts.end(); pi++) {
              Vpt& pt = pts[*pi];
              Point3d pos(pt.x, pt.y, pt.z);
              v3.push_back(pos);
              if (v3.size() == 3)
                break; //we've collected 3 points
            }
            //compute normal
            Vector3d n = ((v3[1] - v3[0]) % (v3[2] - v3[0])).normalize();
            //copy normals
            pi = i.pts.begin();
            for (; pi != i.pts.end(); pi++) {
              Vpt& pt = pts[*pi];
              pt.nx += n[0];
              pt.ny += n[1];
              pt.nz += n[2];
            } //end copying normals
          } //end looping polygons
        } else { // use the information provided
          for (polygon& i : polys) {
            auto ni = i.normals.begin();
            auto pi = i.pts.begin();

            for (; pi != i.pts.end(); pi++, ni++) {
              V& n = normals[*ni];
              Vpt& pt = pts[*pi];
              pt.nx += n.x;
              pt.ny += n.y;
              pt.nz += n.z;
            } //end copying normals

          } //end looping polygons
        }

        //normalize
        for (vector<Vpt>::iterator i = pts.begin(); i != pts.end(); i++) {
          i->normalize();
        }
      }

      vector<Vpt> pts;
      vector<Vpt> colors;
      vector<V> textures;
      vector<V> normals;
      list<polygon> polys;
    };

    class MeshReader {
    public:
      virtual ~MeshReader() {
      }
      bool Read(const std::string& filename) {
        this->m_filename = filename;

        ifstream in(m_filename.c_str());
        if (!in.good()) {
          cerr << "Can't open file " << m_filename << endl;
          return false;
        }
        bool r = Read(in);
        in.close();
        return r;
      }
      const objModel& getModel() const {
        return m_data;
      }
      objModel& getModel() {
        return m_data;
      }

      const std::string& GetTexturePath() const { return this->m_mtl_reader.GetTexturePath(); }
      const std::string& GetMaterialPath() const {return this->m_mtl_path; }

    protected:
      MeshReader() {
      }
      virtual bool Read(istream& in)=0;

      string m_filename;
      objModel m_data;

      string m_mtl_path;

      MTLReader m_mtl_reader;
    };

    class objReader: public MeshReader {
    public:
      objReader() {
      }
      virtual ~objReader() {
      }
    protected:
      bool Read(istream& in) override;
    };

    class offReader: public MeshReader {
    public:
      offReader() {
      }
      virtual ~offReader() {
      }
    protected:
      virtual bool Read(istream& in) override {
        string tmp;
        in >> tmp;
        assert(tmp == "OFF");

        int vsize, fsize, tsize;

        // vertices, faces, unknown
        in >> vsize >> fsize >> tsize;

        // read vertices
        for (auto i = 0; i < vsize; ++i) {
          Vpt pt;
          in >> pt.x >> pt.y >> pt.z;
          m_data.pts.push_back(pt);
        }

        // read faces
        for (auto i = 0; i < fsize; ++i) {
          polygon ply;

          int vs, vid;
          in >> vs;
          assert(vs == 3);
          for (auto j = 0; j < 3; ++j) {
            in >> vid;
            ply.pts.push_back(vid);
          }

          m_data.polys.push_back(ply);
        }

        m_data.compute_v_normal();

        return true;
      }
    };

}}//end spaces

#endif //_OBJ_READER_H_
