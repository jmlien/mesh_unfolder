/*
 * objReader.cpp
 *
 *  Created on: Sep 29, 2016
 *      Author: zxi
 */
#include "objReader.h"

#include "util/Path.h"

istream& operator>>(istream& in, V& v) {
  in >> v.x >> v.y >> v.z;
  return in;
}

ostream& operator<<(ostream& out, const V& v) {
  out << v.x << " " << v.y << " " << v.z;
  return out;
}

bool objReader::Read(istream& in) {

  string tmp;

  //read pts
  while (true) {
    in >> tmp;
    if (tmp == "f")
      break;

    if (tmp == "usemtl") {
      //TODO
    }

    if (tmp == "mtllib") {
      string mtl_path;
      in >> mtl_path;
      cout << "mtllib " << mtl_path << std::endl;
      this->m_mtl_paths.push_back(mtl_path);

      // assuming obj and mtl are in the same dir.
      string full_path = masc::path::DirPath(this->m_filename) + "/" + mtl_path;

      this->m_mtl_reader.Read(full_path);
    }

    if (tmp == "v") {
      Vpt pt;
      in >> pt.x >> pt.y >> pt.z;
      m_data.pts.push_back(pt);
    } else if (tmp == "vn") {
      V pt;
      in >> pt;
      m_data.normals.push_back(pt);
    } else if (tmp == "vt") {
      V pt;
      in >> pt.x >> pt.y;
      m_data.textures.push_back(pt);
    }
    getline(in, tmp);
  }

  //read faces
  polygon poly;
  do {

    in >> tmp;

    if (in.eof())
      break;

    if (isdigit(tmp[0])) { //this defines a vertex

      int pos1 = tmp.find('/');
      int pos2 = tmp.rfind('/');

      int id_v = std::atoi(tmp.substr(0, pos1).c_str()) - 1;
      int id_t = std::atoi(tmp.substr(pos1 + 1, pos2).c_str()) - 1;
      int id_n = std::atoi(tmp.substr(pos2 + 1).c_str()) - 1;

      poly.pts.push_back(id_v);
      poly.normals.push_back(id_n);
      poly.textures.push_back(id_t);

    } else if (tmp == "f") {
      m_data.polys.push_back(poly);
      poly.pts.clear();
      poly.normals.clear();
      poly.textures.clear();
    } else {
      getline(in, tmp);
    }

  } while (!in.eof());

  m_data.polys.push_back(poly);
  m_data.compute_v_normal();

  return true;
}
