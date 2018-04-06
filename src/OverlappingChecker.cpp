/*
 * OverlappingChecker.cpp
 *
 *  Created on: Feb 29, 2016
 *      Author: zxi
 */

#include "OverlappingChecker.h"

#include <float.h>
#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>



#include "unfolder.h"

#include "mathtool/Box.h"

#include "itree/RectKD.h"
#include "itree/MiddleStructure.h"

namespace masc {
namespace unfolding {

///////////////////////////////////////////////////////////////////////////////
// BruteForceChecker
///////////////////////////////////////////////////////////////////////////////

double BruteForceChecker::checkOverlapping(const MESH& mesh,
    const Config& config) {

    //cout<<"BruteForceChecker::checkOverlapping"<<endl;

  int count = 0;
  const int F = this->m_unfolder->getModel()->t_size;

  for (int i = 0; i < F; i++) {
    for (int j = i + 1; j < F; j++) {
      if (this->m_unfolder->checkOverlap(i, j)) {
        count++;
      }
    }
  }

  return count;
}

///////////////////////////////////////////////////////////////////////////////
// ITreeChecker
///////////////////////////////////////////////////////////////////////////////
double ITreeChecker::checkOverlapping(const MESH& mesh, const Config& config) {

    //cout<<"ITreeChecker::checkOverlapping"<<endl;

  const int F = this->m_unfolder->getModel()->t_size;

  typedef Interval<EndPoint> Int;
  typedef RectKD<EndPoint, 2> Rect2D;
  typedef MiddleStructure<Rect2D, 2> MTree;
  typedef vector<Rect2D*>::iterator IT2;

  vector<Rect2D*> rects;

  for (int i = 0; i < F; i++) {

    float max_x = -FLT_MAX, min_x = FLT_MAX, max_y = -FLT_MAX, min_y =
    FLT_MAX;

    for (int p = 0; p < 3; p++) {
      const auto& vi1 = mesh[i][p].second;
      const double a[2] = { vi1[0], vi1[2] };
      max_x = (a[0] > max_x) ? a[0] : max_x;
      min_x = (a[0] < min_x) ? a[0] : min_x;
      max_y = (a[1] > max_y) ? a[1] : max_y;
      min_y = (a[1] < min_y) ? a[1] : min_y;
    }

    //
    Rect2D * p = new Rect2D(i);
    p->setInterval(Int(p, EndPoint(min_x), EndPoint(max_x)), 0);
    p->setInterval(Int(p, EndPoint(min_y), EndPoint(max_y)), 1);

    rects.push_back(p);
  }

  RectangleTree<MTree, 2> Tree;
  if (Tree.Build(rects) == false) {
    //failed, use brute force
    cerr
        << "! WARNING: Build interval tree failed. Resort to brute-force method"
        << endl;
    for (auto rect : rects)
      delete (rect);
    return UINT_MAX;
  }

  uint count = 0;
  for (const auto rect : rects) {
    uint i = rect->getVID();
    Tree.query(rect);
    const Rect2D::Intersect& intersections = rect->getIntersections();

    for (uint j : intersections) {
      if (i > j)
        continue; //this will be checked by (j,i) pair
      if (this->m_unfolder->checkOverlap(i, j)) {
        count++;
      }
    }
  }

  //free mem
  for (auto rect : rects)
    delete (rect);

  return count;
}

///////////////////////////////////////////////////////////////////////////////
// PixelChecker
///////////////////////////////////////////////////////////////////////////////
PixelChecker::PixelChecker(Unfolder *unfolder) :
    OverlappingChecker(unfolder) {

  const int WIDTH = 800;
  const int HEIGHT = 800;

  // render the unfolding onto canvas and compute the area
  this->m_canvas = new cv::Mat(WIDTH, HEIGHT, CV_8UC1, cv::Scalar(0));
}

PixelChecker::~PixelChecker() {
  delete this->m_canvas;
  this->m_canvas = nullptr;
}

double PixelChecker::checkOverlapping(const MESH& mesh, const Config& config) {

//    cout<<"PixelChecker::checkOverlapping"<<endl;

  static int id = 0;

  const int F = mesh.size();

  vector<Vector3d> points(F * 3);

  int vid = 0;
  for (const auto& f : mesh) {
    for (auto i = 0; i < 3; ++i) {
      points[vid++] = f[i].second;
    }
  }

  // creating bounding box of the unfolding
  Box3d box;
  box.setFromPoints(points);

  // computing the scale for rendering
  auto dim = box.getDim();
  auto width = dim[0];
  auto height = dim[2];
  auto scale = (double) this->m_canvas->size[0] / max(width, height);

  if(config.shrink) {
    scale *= config.shrink_factor;
  }

  // shift for rendering
  const Vector3d trans = -box.getMin();

  // is this faster ?
  *m_canvas = cv::Mat::zeros(m_canvas->size[0], m_canvas->size[0], CV_8UC1);

  // render the unfolding onto canvas and compute the area
  const cv::Scalar color(255);
  for (const auto& f : mesh) {
    const auto v0 = (f[0].second + trans) * scale;
    const auto v1 = (f[1].second + trans) * scale;
    const auto v2 = (f[2].second + trans) * scale;

    cv::Point points[3] = { cv::Point(v0[0], v0[2]), cv::Point(v1[0], v1[2]),
        cv::Point(v2[0], v2[2]) };

    cv::fillConvexPoly(*m_canvas, points, 3, color);
  }

  //TODO: for debug only
//  cv::imwrite("cv_overlapping_" + std::to_string(++id) + ".png", canvas);

  const double projected_area = cv::countNonZero(*m_canvas);

  const double expected_area = this->m_unfolder->getModel()->surface_area
      * scale * scale;

// = 1.0 no overlap
// > 1.0 has overlap
  const auto ratio = expected_area / projected_area;

  return ratio;
}

} /* namespace unfolding */
} /* namespace masc */
