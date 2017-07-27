/*
 * ConvexHull2D.h
 *
 *  Created on: Nov 1, 2016
 *      Author: zxi
 */

#ifndef SRC_UTIL_CONVEXHULL2D_H_
#define SRC_UTIL_CONVEXHULL2D_H_

#include <vector>
using namespace std;

#include "mathtool/Vector.h"
using namespace mathtool;

namespace masc {
namespace util {

// Compute the convex for a given points on 2D plane.
// input: points on the XZ plane, y is always 0.
// output: points of the hull boundary on the XZ plane.
void ConvexHull2D(const vector<Vector3d>& points, vector<Vector3d>* hull);

// Compute the area of a 2D convex hull.
double ConvexHull2DArea(const vector<Vector3d>& hull);

}
/* namespace util */
} /* namespace masc */

#endif /* SRC_UTIL_CONVEXHULL2D_H_ */
