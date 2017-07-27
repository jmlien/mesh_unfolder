#include "VRMLWriter.h"

#include <fstream>
using namespace std;

#include "mathtool/Vector.h"
using namespace mathtool;

#include "model.h"

namespace masc {
namespace util {

void writeWrlHeader(int nV, int nT, ofstream& out, const Vector3d& diffuse,
    const Vector3d& specular, const Vector3d& emissive, double ambient,
    double shiniess, double transparency) {
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
  out << "          specularColor " << specular[0] << " " << specular[1] << " "
      << specular[2] << std::endl;
  out << "          emissiveColor " << emissive[0] << " " << emissive[1] << " "
      << emissive[2] << std::endl;
  out << "          shininess " << shiniess << std::endl;
  out << "          transparency " << transparency << std::endl;
  out << "        }" << std::endl;
  out << "      }" << std::endl;
  out << "      geometry IndexedFaceSet {" << std::endl;
  out << "        ccw TRUE" << std::endl;
  out << "        solid TRUE" << std::endl;
  out << "        convex FALSE" << std::endl;
}

// Save a segmented model to VMRL based on face's cluster id
void VRMLWriter::save(const std::string& filename, const model* model) const {
  //TODO not implemented yet
}

// Save all models to VMRL format, each model is a component
void VRMLWriter::save(const std::string& filename,
    const std::vector<model*>& models) const {
  std::cout << "- dumping VMRL format to " << filename << std::endl;

  std::ofstream out(filename);

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

  Vector3d diffuse(mathtool::drand48() * 0.8 + 0.2,
      mathtool::drand48() * 0.8 + 0.2, mathtool::drand48() * 0.8 + 0.2);
  Vector3d specular(0.1, 0.1, 0.1);
  Vector3d emissive(0.0, 0.0, 0.0);
  double ambient = 0.0;
  double shiniess = 0.05;
  double transparency = 0.0;

  for (const model* model_ptr : models) {

    const Vector3d diffuse(mathtool::drand48() * 0.9 + 0.1,
        mathtool::drand48() * 0.9 + 0.1, mathtool::drand48() * 0.9 + 0.1);

    out << std::endl;

    writeWrlHeader(model_ptr->v_size, model_ptr->t_size, out, diffuse, specular,
        emissive, ambient, shiniess, transparency);

    // vertices
    out << "        coord DEF co Coordinate {" << std::endl;
    out << "          point [" << std::endl;
    for (int i = 0; i < model_ptr->v_size; ++i) {
      const auto& v = model_ptr->vertices[i];
      out << "                  " << v.p[0] << " " << v.p[1] << " " << v.p[2]
          << "," << std::endl;
    }

    out << "         ]" << std::endl;
    out << "         }" << std::endl;

    // faces
    out << "        coordIndex [ " << std::endl;
    for (int i = 0; i < model_ptr->t_size; ++i) {
      const auto& f = model_ptr->tris[i];
      out << "            " << f.v[0] << ", " << f.v[1] << ", " << f.v[2]
          << ", -1," << std::endl;
    }

    out << "       ]" << std::endl;
    out << "      }" << std::endl;
    out << "    }" << std::endl;
    out << "  ]" << std::endl;
    out << "}" << std::endl;
  }

  out.close();

}
} // namespace util
} // namespace masc
