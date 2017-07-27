/*
 * Geometry.cpp
 *
 *  Created on: Feb 22, 2016
 *      Author: zxi
 */

#include "Geometry.h"

namespace mathtool {

bool isRightTriangle(const Point3d& p1, const Point3d& p2, const Point3d& p3) {
  Vector3d v1 = Vector3d(p1.get());
  Vector3d v2 = Vector3d(p2.get());
  Vector3d v3 = Vector3d(p3.get());
  return isRightTriangle(v1, v2, v3);
}

bool isRightTriangle(const Vector3d& p1, const Vector3d& p2,
    const Vector3d& p3) {

  const double eps = 1e-3;

  double e1 = (p2 - p1).norm();
  double e2 = (p3 - p2).norm();
  double e3 = (p1 - p3).norm();

  if (fabs(e1 * e1 + e2 * e2 - e3 * e3)/e3/e3 < eps
      || fabs(e1 * e1 + e3 * e3 - e2 * e2)/e2/e2 < eps
      || fabs(e2 * e2 + e3 * e3 - e1 * e1)/e1/e1 < eps)
    return true;

  return false;
}

double triangleArea(const Point3d& p1, const Point3d& p2, const Point3d& p3) {
  Vector3d v1 = Vector3d(p1.get());
  Vector3d v2 = Vector3d(p2.get());
  Vector3d v3 = Vector3d(p3.get());
  return triangleArea(v1, v2, v3);
}

double triangleArea(const Vector3d& p1, const Vector3d& p2,
    const Vector3d& p3) {
  double a = (p2 - p1).norm();
  double b = (p3 - p2).norm();
  double c = (p1 - p3).norm();

  double s = (a + b + c) / 2;
  double area = sqrt(s * (s - a) * (s - b) * (s - c));

  return area;
}
}
