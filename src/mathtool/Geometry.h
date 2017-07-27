/*
 * geometry.h
 *
 *  Created on: Feb 22, 2016
 *      Author: zxi
 */

#ifndef SRC_MATHTOOL_GEOMETRY_H_
#define SRC_MATHTOOL_GEOMETRY_H_

#include "Vector.h"
#include "Point.h"

namespace mathtool {
bool isRightTriangle(const Point3d& p1, const Point3d& p2, const Point3d& p3);
bool isRightTriangle(const Vector3d& p1, const Vector3d& p2,
    const Vector3d& p3);

double triangleArea(const Point3d& p1, const Point3d& p2, const Point3d& p3);
double triangleArea(const Vector3d& p1, const Vector3d& p2, const Vector3d& p3);
}

#endif /* SRC_MATHTOOL_GEOMETRY_H_ */
