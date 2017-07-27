/*
 * OverlappingChecker.h
 *
 *  Created on: Feb 29, 2016
 *      Author: zxi
 */

#ifndef SRC_OVERLAPPINGCHECKER_H_
#define SRC_OVERLAPPINGCHECKER_H_

#include "model.h"
#include "config.h"

// forward declaration
class Unfolder;


namespace cv {
// forward declaration
class Mat;
}

namespace masc {
namespace unfolding {

class OverlappingChecker {
public:
  OverlappingChecker(Unfolder* unfolder) :
      m_unfolder(unfolder) {

  }

  virtual ~OverlappingChecker() {
    m_unfolder = nullptr;
  }

  // check overlapping, greater value means more overlaps
  virtual double checkOverlapping(const MESH& mesh, const Config& config) = 0;
protected:
  // does not own the object
  Unfolder* m_unfolder;
};

// check every pair of facets
class BruteForceChecker: public OverlappingChecker {
public:
  BruteForceChecker(Unfolder* unfolder) :
      OverlappingChecker(unfolder) {

  }
  virtual ~BruteForceChecker() {

  }

  virtual double checkOverlapping(const MESH& mesh, const Config& config)
      override;
};

// use KDTree to speed up
class ITreeChecker: public OverlappingChecker {
public:
  ITreeChecker(Unfolder *unfolder) :
      OverlappingChecker(unfolder) {

  }

  virtual ~ITreeChecker() {

  }

  virtual double checkOverlapping(const MESH& mesh, const Config& config)
      override;
};

// approximate checker by rendering the unfolding and compare the surface area
class PixelChecker: public OverlappingChecker {
public:
  PixelChecker(Unfolder *unfolder);

  virtual ~PixelChecker();

  virtual double checkOverlapping(const MESH& mesh, const Config& config)
      override;

protected:
  cv::Mat* m_canvas;

};

} /* namespace unfolding */
} /* namespace masc */

#endif /* SRC_OVERLAPPINGCHECKER_H_ */
