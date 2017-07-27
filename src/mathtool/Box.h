/*
 * Box.h
 *
 *  Created on: Feb 29, 2016
 *      Author: zxi
 */

#ifndef SRC_MATHTOOL_BOX_H_
#define SRC_MATHTOOL_BOX_H_

#include <vector>
#include <string>

#include "mathtool/Point.h"
#include "mathtool/Vector.h"

namespace mathtool {

class Box2d {
public:
  Box2d();

  Box2d& setFromPoints(const std::vector<Vector2d>& points);

  bool intersect(const Box2d& box) const;

  double x { 0 };
  double y { 0 };
  double width { 0 };
  double height { 0 };
};

class Box3d {
public:
  Box3d();

  Box3d(const std::vector<Point3d>& points) :
      m_empty(false) {
    this->setFromPoints(points);
  }

  Box3d(const std::vector<Vector3d>& points) :
      m_empty(false) {
    this->setFromPoints(points);
  }

  ~Box3d();

  Box3d& addPoint(const Vector3d& p);
  Box3d& addPoint(const Point3d& p);

  Box3d& set(const Vector3d& min, const Vector3d& max) {
    this->m_min = min;
    this->m_max = max;
    return *this;
  }

  Box3d& setFromPoints(const std::vector<Point3d>& points);
  Box3d& setFromPoints(const std::vector<Vector3d>& points);

  const Vector3d& getMin() const {
    return m_min;
  }

  const Vector3d& getMax() const {
    return m_max;
  }

  const Vector3d getDim() const {
    return m_max - m_min;
  }

  const double getSumOfDim() const {
    auto d = this->getDim();
    return d[0] + d[1] + d[2];
  }

private:
  bool m_empty;
  Vector3d m_min;
  Vector3d m_max;
};

} /* namespace mathtool */

#endif /* SRC_MATHTOOL_BOX_H_ */
