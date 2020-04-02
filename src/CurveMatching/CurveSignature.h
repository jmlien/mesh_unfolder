/*
 *  CurveSignature.h
 *  CurveMatching
 *
 *  Created by Roy Shilkrot on 12/7/12.
 *  Copyright (c) 2013 MIT
 *
 *  Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
 *  The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
 *
 *  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 *
 *
 */

#pragma once
#include <opencv2/core/core.hpp>
#include <opencv2/features2d/features2d.hpp>
#include <sys/stat.h>

#include "CurveMatching/opencvstd.h"
using namespace cv;
//


#ifdef ENABLE_PROFILE
#define CV_PROFILE_MSG(msg,code)	\
{\
std::cout << msg << " ";\
double __time_in_ticks = (double)cv::getTickCount();\
{ code }\
std::cout << "DONE " << ((double)cv::getTickCount() - __time_in_ticks)/cv::getTickFrequency() << "s" << std::endl;\
}

#define CV_PROFILE(code)	\
{\
std::cout << #code << " ";\
double __time_in_ticks = (double)cv::getTickCount();\
{ code }\
std::cout << "DONE " << ((double)cv::getTickCount() - __time_in_ticks)/cv::getTickFrequency() << "s" << std::endl;\
}
#else
#define CV_PROFILE_MSG(msg,code) code
#define CV_PROFILE(code) code
#endif


class CurveMatcher
{
public:

	//in this class position type is int for measuing vertex location
	typedef int position_type;

	CurveMatcher()
	{
		RESAMPLE_SIZE=200;           //curve resample size
		SMALLEST_CURVE_SIZE=50;      //smalles curve in the segment DB
		LONGEST_CURVE_SIZE = RESAMPLE_SIZE / 2;
		OFFSET_STEP=2;               //curve segment DB offset incremental size
		CURVATURE_FILTER_SIZE=200;   //number of matches that will pass through curvature filter
	}

	//---
	// set source/target of the match
	//---
	void setSource(const string& filename);
	void setSource(const vector<cv::Point>& a);
	void setSource(const vector<cv::Point2d>& a);
	vector<cv::Point2d>& getSource(){ return source; }

	void setTarget(const string& filename);
	void setTarget(const vector<cv::Point>& b);
	void setTarget(const vector<cv::Point2d>& b);
	vector<cv::Point2d>& getTarget(){ return target; }


	// THE MAIN Function: match the source to the target
	void CompareCurvesUsingSignatureDB(int& a_len, int& a_off, int& b_len, int& b_off, double& score)
	{
		vector<cv::Point2d> ad; ConvertCurve(source, ad);
		vector<cv::Point2d> bd; ConvertCurve(target, bd);
		CompareCurvesUsingSignatureDB(ad, bd, a_len, a_off, b_len, b_off, score);
	}

	void CompareCurvesUsingSignatureDB(
		vector<cv::Point>& a_DB_params, vector<vector<double> > & a_DB,
		vector<cv::Point>& b_DB_params, vector<vector<double> > & b_DB,
		int& a_len,
		int& a_off,
		int& b_len,
		int& b_off,
		double& score
		);


	template<typename T>
	void PrepareSignatureDB(const vector<cv::Point_<T> >& curve, vector<vector<double> >& DB, vector<cv::Point>& DB_params) {
		vector<cv::Point2d> curved; ConvertCurve(curve, curved);
		PrepareSignatureDB(curved, DB, DB_params);
	}


	// prepare curve segment database for comparisons among the shadwos casted by the same model
	// In this case, we would like the min_curve at least half of the given curve.
	template<typename T>
	void PrepareSignatureDB_For_IntraModel_Comparison(const vector<cv::Point_<T> >& curve, vector<vector<double> >& DB, vector<cv::Point>& DB_params)
	{
		vector<cv::Point2d> curved; ConvertCurve(curve, curved);
		PrepareSignatureDB_For_IntraModel_Comparison(curved, DB, DB_params);
	}

	//given a match compute transform
	void computeTransform(int a_len, int a_off, int b_len, int b_off, cv::Mat& transform);
	void computeTransform(vector<cv::Point2d>& a_subset, vector<cv::Point2d>& b_subset, cv::Mat& transform);

	//visualize the match
	void visualizeMatching(int a_len, int a_off, int b_len, int b_off);
	//render matching into an image with name "filename"
	void renderMatching(const string& filename, int a_len, int a_off, int b_len, int b_off);
	void drawMatching(cv::Mat& outout, vector<cv::Point2d> & a_subset, vector<cv::Point2d> & b_subset);

	//compute RMSE of the match
	double ComputeRMSE(int a_len, int a_off, int b_len, int b_off);
	double ComputeRMSE(vector<cv::Point2d>& a_subset, vector<cv::Point2d>& b_subset);

	//compute RMSE of the whole curve using the match
	double ComputeWholeRMSE(int a_len, int a_off, int b_len, int b_off);

	//
	//
	// Access functions!
	//
	//

	//get the length of a given curve
	//in this case, all given curves will be resampled into RESAMPLE_SIZE vertices
	//so the length is always RESAMPLE_SIZE
	int getLength(const vector<cv::Point >& curve) { return RESAMPLE_SIZE; }

	void setCurveResampleSize(int size)
	{
		if (size <= 0){
			cerr << "! Warning: CurveMatcher::setCurveResampleSize has size < 0" << endl; return;
		}
		RESAMPLE_SIZE = size;
	}
	int getCurveResampleSize() const { return RESAMPLE_SIZE; }

	void setSmallestCurveSegmentSize(int size)
	{
		if (size <= 0){
			cerr << "! Warning: CurveMatcher::setSmallestCurveSegmentSize has size < 0" << endl; return;
		}
		SMALLEST_CURVE_SIZE = size;      //smalles curve in the segment DB
	}
	int getSmallestCurveSegmentSize() const { return SMALLEST_CURVE_SIZE; }

	void setLongestCurveSegmentSize(int size)
	{
		if (size <= 0 || size >= RESAMPLE_SIZE)
		{
			cerr << "! Warning: CurveMatcher::setLongestCurveSegmentSize has size smaller than 0 or greater than " << RESAMPLE_SIZE << endl; return;
		}

		LONGEST_CURVE_SIZE = size;
	}
	int getLongestCurveSegmentSize() const { return LONGEST_CURVE_SIZE; }

	void setOffsetStepSize(int size)
	{
		if (size <= 0){
			cerr << "! Warning: CurveMatcher::setOffsetStepSize has size < 0" << endl; return;
		}
		OFFSET_STEP = size;               //curve segment DB offset incremental size
	}

	int getOffsetStepSize() const { return OFFSET_STEP; }

	void setCurvatureFilterSize(int size)
	{
		if (size <= 0) return;
		CURVATURE_FILTER_SIZE = size;   //number of matches that will pass through curvature filter
	}
	int getCurvatureFilterSize() const { return CURVATURE_FILTER_SIZE; }

	//get a subset of curve
	template<typename PT>
	void getSubset(const vector<PT>& curve, int arclen, int offset, vector<PT>& subcurve)
	{
		//make sure that the size of the curve is the same as the RESAMPLE_SIZE
		if (curve.size() != RESAMPLE_SIZE)
		{
			vector<PT> tmp;
			ResampleCurve(curve, tmp, RESAMPLE_SIZE, false);
			subcurve.insert(subcurve.end(), tmp.begin() + offset, tmp.begin() + offset + arclen);
		}
		else
		{
			subcurve.insert(subcurve.end(), curve.begin() + offset, curve.begin() + offset + arclen);
		}
		//
	}

	//load curve from image
	template <typename T>
	void GetCurveForImage(const cv::Mat& filename, vector<cv::Point_<T> >& curve, bool onlyUpper = true, bool getLower = false)
	{
		vector<cv::Point> curve_2i;
		GetCurveForImage(filename, curve_2i, onlyUpper, getLower);
		ConvertCurve(curve_2i, curve);
	}

	//this prepares the curve segment database for both a and b
	void CompareCurvesUsingSignatureDB
		(const vector<cv::Point>& a_DB_params,
		const cv::Point& b_DB_params,
		const vector<vector<double> >& a_DB,
		const vector<double>& b_DB,
		cv::DMatch& best_match);

	void CompareCurvesUsingSignatureDB_curvature_only
		(const vector<cv::Point>& a_DB_params,
		const cv::Point& b_DB_params,
		const vector<vector<double> >& a_DB,
		const vector<double>& b_DB,
		cv::DMatch& best_match);


	//this prepares the curve segment database for both a and b
	void CompareCurvesUsingSignatureDB(const vector<cv::Point>& a_DB_params,
		const vector<cv::Point>& b_DB_params,
		const vector<vector<double> >& a_DB,
		const vector<vector<double> >& b_DB,
		vector<pair<double, cv::DMatch> >& scores_to_matches
		);

	//this prepares the curve segment database for both a and b
	void CompareCurvesUsingSignatureDB(const vector<cv::Point2d>& a,
		const vector<cv::Point2d>& b,
		int& a_len,
		int& a_off,
		int& b_len,
		int& b_off,
		double& score
		);

protected:

	bool fileExists(const std::string& filename);

	template<typename T>
	int closestPointOnCurveToPoint(const vector<cv::Point_<T> >& _tmp, const cv::Point& checkPt, const T cutoff) {
		vector<cv::Point_<T> > tmp = _tmp;
		Mat(tmp) -= Scalar(checkPt.x, checkPt.y);
		vector<float> tmpx, tmpy, tmpmag;
		PolyLineSplit(tmp, tmpx, tmpy);
		magnitude(tmpx, tmpy, tmpmag);
		double minDist = -1;
		cv::Point minLoc; minMaxLoc(tmpmag, &minDist, 0, &minLoc);
		if (minDist < cutoff)
			return minLoc.x;
		else
			return -1;
	}

	template<typename T>
	void saveCurveToFile(const vector<cv::Point_<T> >& curve) {
		static int curve_id = 0;

		stringstream ss; ss << "curves/curves_" << (curve_id++) << ".txt";
		while (fileExists(ss.str())) {
			ss.str("");
			ss << "curves/curves_" << (curve_id++) << ".txt";
		}

		ofstream ofs(ss.str().c_str());
		ofs << curve.size() << "\n";
		for (int i = 0; i < curve.size(); i++) {
			ofs << curve[i].x << " " << curve[i].y << "\n";
		}
		cout << "saved " << ss.str() << "\n";
	}

	template<typename T>
	vector<cv::Point_<T> > loadCurveFromFile(const string& filename) {
		vector<cv::Point_<T> > curve;
		ifstream ifs(filename.c_str());
		int curve_size; ifs >> skipws >> curve_size;
		while (!ifs.eof()) {
			T x, y;
			ifs >> x >> y;
			curve.push_back(cv::Point_<T>(x, y));
		}
		return curve;
	}


	template<typename V>
	cv::Mat_<double> Find2DRigidTransform(const vector<cv::Point_<V> >& a, const vector<cv::Point_<V> >& b,
		cv::Point_<V>* diff = 0, V* angle = 0, V* scale = 0) {
		//use PCA to find relational scale
		Mat_<V> P; Mat(a).reshape(1, a.size()).copyTo(P);
		Mat_<V> Q; Mat(b).reshape(1, b.size()).copyTo(Q);
		PCA a_pca(P, Mat(), PCA::DATA_AS_ROW), b_pca(Q, Mat(), PCA::DATA_AS_ROW);
		double s = sqrt(b_pca.eigenvalues.at<V>(0)) / sqrt(a_pca.eigenvalues.at<V>(0));

		if (s != s)
		{
			for (int i = 0; i < a.size(); i++) if (a[i].x != a[i].x) cout << "a(" << i << ")=" << a[i] << endl;
			//cout << "a_pca=" << a_pca.eigenvalues.at<V>(0) << endl;
		}
		//	cout << a_pca.eigenvectors << endl << a_pca.eigenvalues << endl << a_pca.mean << endl;
		//	cout << b_pca.eigenvectors << endl << b_pca.eigenvalues << endl << b_pca.mean << endl;

		//convert to matrices and subtract mean
		//	Mat_<double> P(a.size(),2),Q(b.size(),2);
		Scalar a_m = Scalar(a_pca.mean.at<V>(0), a_pca.mean.at<V>(1));
		Scalar b_m = Scalar(b_pca.mean.at<V>(0), b_pca.mean.at<V>(1));
		//	for (int i=0; i<a.size(); i++) { P(i,0) = a[i].x - a_m[0]; P(i,1) = a[i].y - a_m[1]; }
		//	for (int i=0; i<b.size(); i++) { Q(i,0) = b[i].x - b_m[0]; Q(i,1) = b[i].y - b_m[1]; }
		P -= repeat((Mat_<V>(1, 2) << a_m[0], a_m[1]), P.rows, 1);
		Q -= repeat((Mat_<V>(1, 2) << b_m[0], b_m[1]), Q.rows, 1);

		//    cout << "new mean for a " << mean(P) << "\n";

		//from http://en.wikipedia.org/wiki/Kabsch_algorithm
		Mat_<double> A = P.t() * Q;
		SVD svd(A);
		Mat_<double> C = svd.vt.t() * svd.u.t();
		double d = (determinant(C) > 0) ? 1 : -1;
		Mat_<double> R = svd.vt.t() * (Mat_<double>(2, 2) << 1, 0, 0, d) * svd.u.t();
		Mat_<double> T = (Mat_<double>(3, 3) << 1, 0, b_m[0] / s, 0, 1, b_m[1] / s, 0, 0, 1) *
			(Mat_<double>(3, 3) << s, 0, 0, 0, s, 0, 0, 0, s) *
			(Mat_<double>(3, 3) << R(0, 0), R(0, 1), 0, R(1, 0), R(1, 1), 0, 0, 0, 1) *
			(Mat_<double>(3, 3) << 1, 0, -a_m[0], 0, 1, -a_m[1], 0, 0, 1)
			;
		if (diff != NULL) {
			diff->x = b_m[0] - a_m[0];
			diff->y = b_m[1] - a_m[1];
		}
		if (angle != NULL) {
			*angle = atan2(R(1, 0), R(0, 0));
		}
		if (scale != NULL) {
			*scale = s;
		}

		return T(Range(0, 2), Range::all());
	}

	template<typename T, typename V>
	cv::Mat_<T> ConvertToMat(const vector<vector<V> >& mnt_DB) {
		Mat_<T> mnt_DB_m(mnt_DB.size(), mnt_DB[0].size());
		for (int i = 0; i < mnt_DB.size(); i++) {
			for (int j = 0; j < mnt_DB[i].size(); j++) {
				mnt_DB_m(i, j) = (T)(mnt_DB[i][j]);
			}
		}
		return mnt_DB_m;
	}

	template<typename T>
	void drawCurvePoints(cv::Mat& img, const vector<cv::Point_<T> >& curve_, const cv::Scalar& color, int thickness) {
		vector<cv::Point> curve;
		ConvertCurve(curve_, curve);
		for (int i = 0; i < curve.size(); i++) {
			circle(img, curve[i], 3, color, thickness);
		}
	}

	void GetCurveForImage(const cv::Mat& filename, vector<cv::Point>& curve, bool onlyUpper = true, bool getLower = false);
	void GetCurveForImage(const cv::Mat& filename, vector<cv::Point>& whole, vector<cv::Point>& curve_upper, vector<cv::Point>& curve_lower);

	template<int x, int y>
	void imshow_(const std::string& str, const cv::Mat& img)
	{
		cv::Mat big;
		resize(img, big, cv::Size(x, y), -1, -1, INTER_NEAREST);
		imshow(str, big);
	}

	template<typename T, typename V>
	vector<cv::Point_<V> > YnormalizedCurve(const vector<cv::Point_<T> >& curve) {
		vector<T> curvex, curvey; PolyLineSplit(curve, curvex, curvey);
		double minval, maxval;
		minMaxIdx(curvey, &minval, &maxval);
		vector<cv::Point_<V> > curveout;
		for (int i = 0; i < curve.size(); i++) {
			curveout.push_back(cv::Point_<V>(i, (curvey[i] - minval) / (maxval - minval)));
		}
		return curveout;
	}

	cv::Mat_<double> GetSmithWatermanHMatrix(const vector<pair<char, int> >& a, const vector<pair<char, int> >& b);

	double MatchSmithWaterman(const vector<pair<char, int> >& a, const vector<pair<char, int> >& b, vector<cv::Point>& matching);

	// prepare curve segment database for generic purposes
	// in particular for comparison between different models
	void PrepareSignatureDB(const vector<cv::Point2d>& curve, vector<vector<double> >& DB, vector<cv::Point>& DB_params);


	// prepare curve segment database for comparisons among the shadwos casted by the same model
	// In this case, we would like the min_curve at least half of the given curve.
	void PrepareSignatureDB_For_IntraModel_Comparison(const vector<cv::Point2d>& curve, vector<vector<double> >& DB, vector<cv::Point>& DB_params);


private:

	vector<cv::Point2d> source;
	vector<cv::Point2d> target;

	//parameters used in this matcher
	int RESAMPLE_SIZE;            //curve resample size
	int SMALLEST_CURVE_SIZE;      //smallest curve in the segment DB
	int LONGEST_CURVE_SIZE;       //longest curve in the segment DB
	int OFFSET_STEP;              //curve segment DB offset incremental size
	int CURVATURE_FILTER_SIZE;    //number of matches that will pass through curvature filter
};

//a dummy function for compatiability
inline float getArcLength(const vector<cv::Point>& curve)
{
	return (float)curve.size();
}
