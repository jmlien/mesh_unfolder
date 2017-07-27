#ifndef _UTIL_VMRL_WRITER_H_
#define _UTIL_VMRL_WRITER_H_

#include <iostream>
#include <vector>
#include <string>

struct model;

namespace masc {
namespace util {

class VRMLWriter {
public:
  VRMLWriter() {
  }
  ~VRMLWriter() {
  }

  // Save a segmented model to VMRL based on face's cluster id
  void save(const std::string& filename, const model* model) const;

  // Save all models to VMRL format, each model is a component
  void save(const std::string& filename,
      const std::vector<model*>& models) const;
};

} //namespace util
} //namespace masc

#endif //_UTIL_VMRL_WRITER_H_
