/*
 * Box.cpp
 *
 *  Created on: Feb 29, 2016
 *      Author: zxi
 */

#include "Box.h"

#include <float.h>

namespace mathtool {

Box2d::Box2d() {

}

Box2d& Box2d::setFromPoints(const std::vector<Vector2d>& points) {
  Vector2d min(FLT_MAX, FLT_MAX);
  Vector2d max(FLT_MIN, FLT_MIN);

  for (const auto& p : points) {
    for (int i = 0; i < 2; ++i) {
      if (p[i] < min[i])
        min[i] = p[i];
      if (p[i] > max[i])
        max[i] = p[i];
    }
  }

  this->x = min[0];
  this->y = min[1];
  this->width = max[0] - min[0];
  this->height = max[1] - min[1];

  return *this;
}

bool Box2d::intersect(const Box2d& box) const {
  return (fabs(x - box.x) * 2 < (this->width + box.width))
      && (fabs(y - box.y) * 2 < (this->height + box.height));
}

Box3d::Box3d() {
  m_empty = true;
}

Box3d::~Box3d() {
}

Box3d& Box3d::addPoint(const Vector3d& p) {
  if (this->m_empty) {
    this->m_min = this->m_max = p;
    this->m_empty = false;
  } else {
    for (int i = 0; i < 3; ++i) {
      if (p[i] < m_min[i])
        m_min[i] = p[i];
      if (p[i] > m_max[i])
        m_max[i] = p[i];
    }
  }
  return *this;

}
Box3d& Box3d::addPoint(const Point3d& p) {
  if (this->m_empty) {
    this->m_min = this->m_max = Vector3d(p.get());
    this->m_empty = false;
  } else {
    for (int i = 0; i < 3; ++i) {
      if (p[i] < m_min[i])
        m_min[i] = p[i];
      if (p[i] > m_max[i])
        m_max[i] = p[i];
    }
  }
  return *this;
}

Box3d& Box3d::setFromPoints(const std::vector<Point3d>& points) {
  Vector3d min(FLT_MAX, FLT_MAX, FLT_MAX);
  Vector3d max(FLT_MIN, FLT_MIN, FLT_MIN);

  for (const auto& p : points) {
    for (int i = 0; i < 3; ++i) {
      if (p[i] < min[i])
        min[i] = p[i];
      if (p[i] > max[i])
        max[i] = p[i];
    }
  }

  return this->set(min, max);
}

Box3d& Box3d::setFromPoints(const std::vector<Vector3d>& points) {
  Vector3d min(FLT_MAX, FLT_MAX, FLT_MAX);
  Vector3d max(FLT_MIN, FLT_MIN, FLT_MIN);

  for (const auto& p : points) {
    for (int i = 0; i < 3; ++i) {
      if (p[i] < min[i])
        min[i] = p[i];
      if (p[i] > max[i])
        max[i] = p[i];
    }
  }

  return this->set(min, max);
}

} /* namespace mathtool */
