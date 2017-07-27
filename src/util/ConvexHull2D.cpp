/*
 * ConvexHull2D.cpp
 *
 *  Created on: Nov 1, 2016
 *      Author: zxi
 */

#include "ConvexHull2D.h"

#include <opencv2/imgproc/imgproc.hpp>

namespace masc {
namespace util {

void ConvexHull2D(const vector<Vector3d>& points, vector<Vector3d>* hull) {

  hull->clear();

  vector<cv::Point2f> input;
  vector<int> hull_indices;

  for (const auto& pt : points)
    input.push_back(cv::Point2f(pt[0], pt[2]));

  cv::convexHull(input, hull_indices);

  for (const int index : hull_indices) {
    hull->push_back(points[index]);
  }
}

double ConvexHull2DArea(const vector<Vector3d>& hull) {
  double area = 0.0;
  for (int i = 0; i < hull.size(); ++i) {
    int j = (i + 1) % hull.size();
    area += (hull[i][0] * hull[j][2] - hull[j][0] * hull[i][2]);
  }
  return area * 0.5;
}

} /* namespace util */
} /* namespace masc */
