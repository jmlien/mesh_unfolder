/*
 * Path.h
 *
 *  Created on: Sep 29, 2016
 *      Author: zxi
 */

#ifndef SRC_UTIL_PATH_H_
#define SRC_UTIL_PATH_H_

#include <string>

namespace masc {
namespace path {

// Get the directory path
// /path/to/file.txt -> /path/to
std::string DirPath(const std::string& path);

// Get the base name of a path
// /path/to/file.txt -> file.txt
std::string Basename(const std::string& path);

// Get the stem of a path
// /path/to/file.txt -> file
std::string Stem(const std::string& path);

// Get the stem of a path
// /path/to/file.txt -> .txt
// /path/to/file -> "" empty string
std::string Extension(const std::string& path);

// Check whether a path exists or not.
bool Exists(const std::string& path);

} /* namespace path */
} /* namespace masc */

#endif /* SRC_UTIL_PATH_H_ */
