//create a list of shadows for each model

#pragma once

#include <map>

//openCV
#include "opencv2/core/core.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include "CurveMatching/curve_db_param.h"

//own stuff
#include "polygon/polygon.h"
#include "CurveMatching/img2ply.h"
#include "CurveMatching/opencvstd.h"
#include "CurveMatching/CurveCSS.h"
#include "CurveMatching/CurveSignature.h"
#include "CurveMatching/ArcLengthCurveSignature.h"
//#include "simple_svg_export.h"
///-------

//
#include "Vector.h"

//-----------------------------------------------------------------------------
//
//
// debug parameters in this file
//
//
//-----------------------------------------------------------------------------

#define SHOW_DIST_TRANSFORM         00

// Initial directions for creating shadows of a given model, one direction per shadow.
// these shadows will later be clustered to create representative shadow
//

//-----------------------------------------------------------------------------

//convert from c_ply to contour date strcuture (a vector of points)
template<class PT> void ply2contour(const c_ply& ply, vector<PT>& contour) {
	//contour.reserve(ply.getSize());
	ply_vertex * ptr = ply.getHead();
	do {
		const mathtool::Point2d& pt = ptr->getPos();
		contour.push_back(PT(pt[0], pt[1]));
		ptr = ptr->getNext();
	} while (ptr != ply.getHead());
}

//-------------------------------------------------------------------------------------------------------------------------------
template<typename MATCHER, typename PT>
struct CSShape {

	//matching database in opencv format...
	void buildCurveSegmentDB()
	{
		MATCHER matcher;
		CURVE_DB_PARAM<MATCHER>::initializeMatcher(matcher);

		int csize = contours.size();

		for (int i = 0; i < csize; i++) {
			for (c_polygon::iterator k = contours[i].begin();
					k != contours[i].end(); k++) {
				contourDBs.push_back(contour_DB());
				contour_DB& db = contourDBs.back();
				ply2contour(*k, db.cv_contour);

				//what is this??
				db.cv_contour.push_back(db.cv_contour.front());  //first element
                //db.cv_contour.push_back(*(++db.cv_contour.begin())); //second element...
				//

				matcher.PrepareSignatureDB(db.cv_contour, db.curvatures, db.curve_segments);

				//cout << "length=" << matcher.getLength(db.cv_contour) << endl;

				//db.curve_segments.back().x is the length of the longest segment created in the DB
				//db.computeConflictList(matcher.getLength(db.cv_contour)); // db.curve_segments.back().x);
			} //end for k (each contour, including external and hole boundaries)

		}				//end for i (each contour polygon)

		//cout << "There are " << contourDBs.size() << " contourDBs" << endl;
	}

	CSShape& operator=(CSShape& other)
	{
		//in contourDBs, each element (contour_DB) is created for each c_polygon in contours
		this->contourDBs = other.contourDBs; //note: the size of contourDBs may be different from contours size (why?)
		domain=other.domain.clone();
		distfield=other.distfield.clone(); //distance field of the image above
		contours=other.contours;
		return *this;
	}

	//
	//check if a given c_ply with its transform is inside this target shadow or not
	//there are two parameters used to control how relaxed the inside/outside definition should be
	//tolerable_max_diff: controls the max distance difference between the target and shadow boundary
	//tolerable_sum_diff: controls the total, i.e. area, that the shadow is outside the target
	//
	bool inside(const vector<PT>& ply, cv::Mat& tranform,
			float tolerable_max_diff = 10, float tolerable_sum_diff = FLT_MAX);

	//
	//compute distance transform on target_image using open cv function
	//
	void builDistField()
	{
		cv::distanceTransform(this->domain, distfield, DIST_L2, DIST_MASK_5);

#if SHOW_DIST_TRANSFORM
		{
			double min;
			double max;
			cv::minMaxIdx(distfield, &min, &max);
			cv::Mat adjMap;
			cv::convertScaleAbs(distfield, adjMap, 255 / max);
			cv::applyColorMap(adjMap, adjMap, cv::COLORMAP_JET);
			cv::imshow("distance field", adjMap);
			cv::waitKey();
		}
#endif
	}

	void draw(cv::Mat& mat, cv::Scalar& color) {
		vector<vector<cv::Point> > cv_contours;
		for (int i = 0; i < contourDBs.size(); i++) {
			auto& target_db = contourDBs[i];
			{
				//TODO: WHY DO WE CONVERT Point2d (target_db.cv_contour) to Point
				//      JML: If it is not Point but Point2d, OpenCV crashes
				vector<cv::Point> tmp2;
				ConvertCurve(target_db.cv_contour, tmp2);
				cv_contours.push_back(tmp2);
			}
		}

		drawContours(mat, cv_contours, -1, color, FILLED);
	}

	struct contour_DB {
		vector<PT> cv_contour; //a contour, mostly type PT is a open_cv format

		//data base signature
		vector<vector<double> > curvatures; //a list of curvatures for each curve segment
		vector<PT> curve_segments; //a list of curve segments, first values (PT.x) is the length, the second value (PT.y) is the offset of the head of polyline

		//overlapping curve segments
		vector<list<unsigned int > > conflict_list; // list<ids> are segments overlapping with each other
	};

	//in contourDBs, each element (contour_DB) is created for each c_polygon in contours
	vector<contour_DB> contourDBs; //note: the size of contourDBs may be different from contours size (why?)

	//cv::Mat image;
	cv::Mat domain;
	cv::Mat distfield; //distance field of the image above
	vector<c_polygon> contours;
};

//check if a given shadow with its transform is inside this shadow or not
template<typename MATCHER, typename PT> bool CSShape<MATCHER, PT>::inside(
		const vector<PT>& ply, cv::Mat& tranform,
		float tolerable_max_diff, float tolerable_sum_diff) {
	//cout << "inside transform: " << endl;
	//cout << tranform ;

	//build a new polygon from shadow.contour;
	//c_polygon trans_contour;
	//transformPolygon(shadow.contour, tranform, trans_contour);

	//JML: Again this is another case that we will have to draw cv::Point instead of cv::Point2d
	//otherwise openCV crashed
	vector<vector<cv::Point> > shadow_trans;
	{
		vector<cv::Point2d> tmp;
		cv::transform(ply, tmp, tranform);
		vector<cv::Point> tmp2;
		ConvertCurve(tmp, tmp2);
		shadow_trans.push_back(tmp2);
	}

	//convert the contour to image
	cv::Mat shadow_img(distfield.size(), CV_8UC1, cv::Scalar::all(255)); //a white image

	//fill in with black color
	drawContours(shadow_img, shadow_trans, 0, cv::Scalar::all(0), FILLED);

	//check if shadow_trans is inside target_contours
	int bad_pixel_count = 0;
	float sum_bad_distance = 0; //sum of distance from the target boundary
	float max_bad_distance = 0; //max of distance from the target boundary

	for (int w = 0; w < distfield.size().width; w++) {
		for (int h = 0; h < distfield.size().height; h++) {

			float d = distfield.at<float>(w, h);
			unsigned char s = shadow_img.at<unsigned char>(w, h);

			if (d > 0 && s == 255)
				continue; //outside target and outside shadow
			if (d == 0 && s == 0)
				continue; //inside target and inside shadow
			if (d == 0 && s == 255)
				continue; //inside target and outside shadow (this is area not covered by the shadow, which is fine)

			//now we have a pixel that is inside the shadow but outside the target...
			bad_pixel_count++;
			sum_bad_distance += d;
			if (d > max_bad_distance)
				max_bad_distance = d;
		} //end h
	} //end w

	//cout << "outside max_bad_distance=" << max_bad_distance << " bad_pixel_count=" << bad_pixel_count << endl;

#if SHOW_DIST_TRANSFORM
	//{//draw distance field
	//	double min;
	//	double max;
	//	cv::minMaxIdx(distfield, &min, &max);
	//	cv::Mat adjMap;
	//	cv::convertScaleAbs(distfield, adjMap, 255 / max);
	//	cv::applyColorMap(adjMap, adjMap, cv::COLORMAP_JET);

	//	this->draw(adjMap, cv::Scalar::all(0));
	//	drawContours(adjMap, shadow_trans, 0, cv::Scalar::all(255), CV_FILLED);
	//	cv::imshow("distance field", adjMap);
	//	cv::waitKey();
	//}
	//cv::imshow("distfield", distfield);
	//cv::imshow("shadow_img", shadow_img);
#endif

	if (max_bad_distance <= tolerable_max_diff
			&& bad_pixel_count <= tolerable_sum_diff) {
		//cout << "inside" << endl;
		return true;
	}

	return false;
}
