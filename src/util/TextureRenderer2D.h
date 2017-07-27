/*
 * TextureRenderer2D.h
 *
 *  Created on: Sep 30, 2016
 *      Author: zxi
 */

#ifndef SRC_UTIL_TEXTURERENDERER2D_H_
#define SRC_UTIL_TEXTURERENDERER2D_H_

#include <string>
#include <opencv2/core/core.hpp>
#include "mathtool/Vector.h"

using namespace mathtool;

namespace masc {
namespace util {

class TextureRenderer2D {
public:
  // Create a TextureRenderer2D with a texture file
  TextureRenderer2D(const std::string& texture_path);
  virtual ~TextureRenderer2D();

  void SetCanvasSize(int width, int height);

  // Given a triangle {p0, p1, p2} on 2D canvas and its UV coordinates {uv0, uv1, uv2},
  // render the texture of the triangle into 2D bitmap whose size equals to its bounding box's.
  bool Render(const Vector2d& p0, const Vector2d& p1, const Vector2d& p2,
      const Vector2d& uv0, const Vector2d& uv1, const Vector2d& uv2);

  // Export the rendering, encoded bitmap as PMG format in base_64 encoding.
  void ExportToPngBase64(string* png_base64);

protected:
  std::string m_texture_path;
  cv::Mat m_texture;
  cv::Mat m_canvas;
  int m_texture_width;
  int m_texture_height;
  int m_id;
  int m_canvas_width;
  int m_canvas_height;
};

} /* namespace util */
} /* namespace masc */

#endif /* SRC_UTIL_TEXTURERENDERER2D_H_ */
