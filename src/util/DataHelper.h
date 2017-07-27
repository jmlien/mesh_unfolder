/*
 * DataHelper.h
 *
 *  Created on: Mar 4, 2016
 *      Author: zxi
 */

#ifndef SRC_UTIL_DATAHELPER_H_
#define SRC_UTIL_DATAHELPER_H_

#include <vector>
#include <fstream>
#include <sstream>

namespace masc {
namespace util {

// read the into a list of given type
template<typename T>
std::vector<T> readList(const std::string& filename) {
  std::vector<T> result;

  std::ifstream fin(filename);
  if (!fin.good()) {
    std::cerr << "!Error! Failed to open " << filename << std::endl;
    return result;
  }

  // Read a line into a local string.
  std::string line;

  while (std::getline(fin, line)) {
    T item;
    // convert the line into a stream and read the key/value from the line
    std::stringstream linestream(line);
    linestream >> item;

    result.push_back(item);

    // If the linestream is bad, then reading the key/value failed
    // If reading one more char from the linestream works then there is extra crap in the line
    // thus we have bad data on a line.
    //
    // In either case set the bad bit for the input stream.
    char c;
    if ((!linestream) || (linestream >> c)) {
      fin.setstate(std::ios::badbit);
    }
  }

  fin.close();

  return result;
}

template<typename T>
void writeList(const std::string& filename, const std::vector<T>& list,
    char separator = '\n') {
  std::ofstream fout(filename);
  if (!fout.good()) {
    std::cerr << "!Error! Failed to open " << filename << std::endl;
    return;
  }

  for (const auto& item : list)
    fout << item << separator;

  fout.close();
}

}
}

#endif /* SRC_UTIL_DATAHELPER_H_ */
