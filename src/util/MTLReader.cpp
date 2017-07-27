/*
 * MTLReader.cpp
 *
 *  Created on: Sep 29, 2016
 *      Author: zxi
 */

#include "MTLReader.h"

#include <fstream>

#include "util/Path.h"

bool MTLReader::Read(const std::string& path) {
  this->m_filename = path;

  std::cout << "MTLReader reading from: " << path << std::endl;

  std::ifstream in(this->m_filename);
  if (!in.good())
    return false;

  return this->Read(in);
}

bool MTLReader::Read(std::istream& in) {
  std::string tmp;
  while (in >> tmp) {
    if (tmp == "map_Ka") {
      in >> this->m_texture_filename;
      this->m_texture_filename = masc::path::DirPath(this->m_filename) + "/"
          + this->m_texture_filename;
      std::cout << "map_Ka " << this->m_texture_filename << std::endl;
    }
  }
  return false;
}
