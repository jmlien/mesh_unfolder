/*
 * TextureRenderer2D.cpp
 *
 *  Created on: Sep 30, 2016
 *      Author: zxi
 */

#include "TextureRenderer2D.h"

#include <iostream>
#include <string>
#include <opencv2/highgui/highgui.hpp>
#include "opencv2/imgproc/imgproc.hpp"
#include "util/Base64.h"
#include "util/Path.h"

using namespace masc::path;

namespace masc {
namespace util {

TextureRenderer2D::TextureRenderer2D(const std::string& texture_path) :
    m_texture_path(texture_path), m_texture_width(0), m_texture_height(0), m_id(
        0), m_canvas_width(0), m_canvas_height(0) {

  std::cerr << "- [TextureRenderer2D::TextureRenderer2D] Reading texture from "
      << texture_path << std::endl;
  m_texture = cv::imread(m_texture_path, CV_LOAD_IMAGE_COLOR);

  m_texture_width = m_texture.size[0];
  m_texture_height = m_texture.size[1];

  std::cerr << "- [TextureRenderer2D::TextureRenderer2D] Texture loaded size: "
      << m_texture_width << "*" << m_texture_height << std::endl;

  assert(m_texture_width > 0 && m_texture_height > 0);

  cv::cvtColor(m_texture, m_texture, CV_BGR2BGRA);
}

TextureRenderer2D::~TextureRenderer2D() {
}

void TextureRenderer2D::SetCanvasSize(int width, int height) {
  m_canvas_width = width;
  m_canvas_height = height;

  // create a transparent canvas...
  m_canvas = cv::Mat(height, width, CV_8UC4, cv::Scalar(0, 0, 0, 0));

  std::cerr << "- [TextureRenderer2D::SetCanvasSize] Canvas size = " << width
      << "*" << height << std::endl;

  // Clean the canvas before drawing
  m_canvas = cv::Scalar(0, 0, 0, 0);
}

bool TextureRenderer2D::Render(const Vector2d& p0, const Vector2d& p1,
    const Vector2d& p2, const Vector2d& uv0, const Vector2d& uv1,
    const Vector2d& uv2) {

  vector<cv::Point2f> srcTri(3);
  vector<cv::Point2f> dstTri(3);

  // Create an temp canvas for drawing
  cv::Mat canvas = cv::Mat(m_canvas_height, m_canvas_width, CV_8UC4,
      cv::Scalar(0, 0, 0, 0));

  /// Set your 3 points to calculate the  Affine Transform
  srcTri[0] = cv::Point2f(uv0[0] * m_texture_width, uv0[1] * m_texture_height);
  srcTri[1] = cv::Point2f(uv1[0] * m_texture_width, uv1[1] * m_texture_height);
  srcTri[2] = cv::Point2f(uv2[0] * m_texture_width, uv2[1] * m_texture_height);

  dstTri[0] = cv::Point2f(p0[0], p0[1]);
  dstTri[1] = cv::Point2f(p1[0], p1[1]);
  dstTri[2] = cv::Point2f(p2[0], p2[1]);

  // coordinates in int
  vector<cv::Point2i> srcTriCropped;

  for (const auto pt : srcTri)
    srcTriCropped.push_back(cv::Point2i(pt.x + 0.5, pt.y + 0.5));

  // Get mask by filling triangle
  cv::Mat mask = cv::Mat(m_texture_height, m_texture_width, CV_8U,
      cv::Scalar(0));
  cv::fillConvexPoly(mask, srcTriCropped, cv::Scalar(255), cv::LINE_AA, 0);

  cv::Mat masked_texture;
  m_texture.copyTo(masked_texture, mask);

  /// Get the Affine Transform
  auto warp_mat = cv::getAffineTransform(srcTri, dstTri);

  /// Apply the Affine Transform just found to the src image
  cv::warpAffine(masked_texture, canvas, warp_mat, canvas.size(),
      cv::INTER_LINEAR, cv::BORDER_TRANSPARENT);

  vector<cv::Point2i> dstTriCropped;
  for (const auto pt : dstTri)
    dstTriCropped.push_back(cv::Point2i(pt.x + 0.5, pt.y + 0.5));

  mask = cv::Mat(m_canvas_height, m_canvas_width, CV_8U, cv::Scalar(0));
  cv::fillConvexPoly(mask, dstTriCropped, cv::Scalar(255), cv::LINE_AA, 0);
  canvas.copyTo(m_canvas, mask);

//  const string dump_path = DirPath(m_texture_path) + "/dump_"
//      + std::to_string(++m_id) + ".png";
//
//  cerr << " - [TextureRenderer2D::Render] Dump image to " << dump_path << endl;
//  cv::imwrite(dump_path, m_canvas);

  return false;
}

void TextureRenderer2D::ExportToPngBase64(string* png_base64) {
  std::vector<uchar> buffer;

  cv::imencode(".png", m_canvas, buffer);

  const auto base64_string = base64_encode(&(*buffer.begin()), buffer.size());

  *png_base64 = "data:image/png;base64," + base64_string;

}

} /* namespace util */
} /* namespace masc */
