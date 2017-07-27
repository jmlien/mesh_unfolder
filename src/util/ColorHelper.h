/*
 * ColorUtil.h
 *
 *  Created on: Mar 6, 2016
 *      Author: zxi
 */

#ifndef SRC_UTIL_COLORHELPER_H_
#define SRC_UTIL_COLORHELPER_H_

#include "mathtool/Vector.h"
using namespace mathtool;

namespace masc {
namespace util {

/*
 Return a RGB colour value given a scalar v in the range [vmin,vmax]
 In this case each colour component ranges from 0 (no contribution) to
 1 (fully saturated), modifications for other ranges is trivial.
 The colour is clipped at the end of the scales if v is outside
 the range [vmin,vmax]
 */

Vector3d getColour(double v, double vmin, double vmax) {
  Vector3d c(1.0, 1.0, 1.0); // white
  double dv;

  if (v < vmin)
    v = vmin;
  if (v > vmax)
    v = vmax;
  dv = vmax - vmin;

  if (v < (vmin + 0.25 * dv)) {
    c[0] = 0;
    c[1] = 4 * (v - vmin) / dv;
  } else if (v < (vmin + 0.5 * dv)) {
    c[0] = 0;
    c[2] = 1 + 4 * (vmin + 0.25 * dv - v) / dv;
  } else if (v < (vmin + 0.75 * dv)) {
    c[0] = 4 * (v - vmin - 0.5 * dv) / dv;
    c[2] = 0;
  } else {
    c[1] = 1 + 4 * (vmin + 0.75 * dv - v) / dv;
    c[2] = 0;
  }

  return c;
}
}
}

#endif /* SRC_UTIL_COLORHELPER_H_ */
