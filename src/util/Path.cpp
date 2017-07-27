/*
 * Path.cpp
 *
 *  Created on: Sep 29, 2016
 *      Author: zxi
 */

#include "Path.h"

#include <fstream>

namespace masc {
namespace path {

std::string DirPath(const std::string& path) {
  size_t pos = path.find_last_of("\\/");
  return (std::string::npos == pos) ? "" : path.substr(0, pos);
}

std::string Basename(const std::string& path) {
  size_t pos = path.find_last_of("\\/");
  return (std::string::npos == pos) ? path : path.substr(pos + 1);
}

bool Exists(const std::string& path) {
  return std::ifstream(path).good();
}

std::string Stem(const std::string& path) {
  std::string base_name = Basename(path);
  auto pos = base_name.find_last_of(".");
  if (pos == std::string::npos)
    return base_name;

  return base_name.substr(0, pos);
}

std::string Extension(const std::string& path) {
  auto pos = path.find_last_of(".");
  if (pos == std::string::npos)
    return "";

  return path.substr(pos);
}

} /* namespace path */
} /* namespace masc */
