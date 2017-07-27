#pragma once

#include "polygon/polygon.h"
#include "CurveMatching/opencvstd.h"

using namespace masc::polygon;

//given a black/whit image with size (widthxheight) return a ploygon representing the contour
//NOTE: this function only extracts the largest boundary and its internal holes
bool img2ply(unsigned char * img, unsigned int width, unsigned int height, c_polygon& polygon);

//given a black/whit image with size (widthxheight) return a ploygon representing the contour
bool img2ply(cv::Mat src, c_polygon& polygon);



//Note: this function extracts all contours except those with size smaller than threshold
bool img2ply_complete(unsigned char * img, unsigned int width, unsigned int height, vector<c_polygon>& polygons, float threshold);


//Note: this function extracts all contours except those with size smaller than threshold
bool img2ply_complete(cv::Mat& img, vector<c_polygon>& polygons, float threshold);
