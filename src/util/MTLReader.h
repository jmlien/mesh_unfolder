/*
 * MTLReader.h
 *
 *  Created on: Sep 29, 2016
 *      Author: zxi
 */

#ifndef SRC_UTIL_MTLREADER_H_
#define SRC_UTIL_MTLREADER_H_

#include <string>
#include <iostream>

class MTLReader {
public:
  MTLReader() {
  }
  virtual ~MTLReader() {
  }
  bool Read(const std::string& path);
  bool Read(std::istream& in);

  const std::string& GetTexturePath() const {
    return m_texture_filename;
  }
protected:
  // full path of the mtl file
  std::string m_filename;
  // full path of the texture file, only supports one.
  std::string m_texture_filename;
};

#endif /* SRC_UTIL_MTLREADER_H_ */
