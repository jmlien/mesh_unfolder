#ifndef ARC_LENGTH_CURVE_SIGNATURE_H_
#define ARC_LENGTH_CURVE_SIGNATURE_H_

#include "opencv2/features2d/features2d.hpp"
#include <sys/stat.h>

#include "CurveMatching/opencvstd.h"
using namespace cv;

#include "CurveMatching/CurveCSS.h"

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

float getArcLength(const vector<cv::Point2d>& curve);

//
// this is arc length based matcher
// all matches are between segments with the same length only
//
class CurveArcLengthMatcher
{
public:

	//in this class position type is float for measuing arc length
	typedef float position_type;

	CurveArcLengthMatcher()
	{
		RESAMPLE_SIZE=200;           //curve resample size
		SMALLEST_CURVE_SIZE=50;      //smallest curve in the segment DB
		LONGEST_CURVE_SIZE = 100;    //longest curve in the segment DB measured in arc length
		OFFSET_STEP=2;               //curve segment DB offset incremental size
		CURVATURE_FILTER_SIZE=200;   //number of matches that will pass through curvature filter
	}

	//---
	// set source/target of the match
	//---
	void setSource(const string& filename);
	void setSource(const vector<cv::Point2d>& a);
	vector<cv::Point2d>& getSource(){ return source; }

	void setTarget(const string& filename);
	void setTarget(const vector<cv::Point2d>& b);
	vector<cv::Point2d>& getTarget(){ return target; }

	//-------------------------------------------------------------------------

	// THE MAIN Function: match the source to the target
	void CompareCurvesUsingSignatureDB(float& a_len, float& a_off, float& b_len, float& b_off, double& score)
	{
		CompareCurvesUsingSignatureDB(source, target, a_len, a_off, b_len, b_off, score);
	}

	//this prepares the curve segment database for both a and b
	void CompareCurvesUsingSignatureDB
		(const vector<cv::Point2d>& a_DB_params,
		const cv::Point2d& b_DB_params,
		const vector<vector<double> >& a_DB,
		const vector<double>& b_DB,
		cv::DMatch& best_match);

	void CompareCurvesUsingSignatureDB_curvature_only
		(const vector<cv::Point2d>& a_DB_params,
		const cv::Point2d& b_DB_params,
		const vector<vector<double> >& a_DB,
		const vector<double>& b_DB,
		cv::DMatch& best_match);

	//this prepares the curve segment database for both a and b
	void CompareCurvesUsingSignatureDB(const vector<cv::Point2d>& a_DB_params,
		const vector<cv::Point2d>& b_DB_params,
		const vector<vector<double> >& a_DB,
		const vector<vector<double> >& b_DB,
		vector<pair<double, cv::DMatch> >& scores_to_matches
		);

	//this prepares the curve segment database for both a and b
	void CompareCurvesUsingSignatureDB(const vector<cv::Point2d>& a,
		const vector<cv::Point2d>& b,
		float& a_len,
		float& a_off,
		float& b_len,
		float& b_off,
		double& score
		);

	void CompareCurvesUsingSignatureDB(
		vector<cv::Point2d>& a_DB_params, vector<vector<double> > & a_DB,
		vector<cv::Point2d>& b_DB_params, vector<vector<double> > & b_DB,
		float& a_len,
		float& a_off,
		float& b_len,
		float& b_off,
		double& score
		);


	//-------------------------------------------------------------------------

	template<typename T>
	void PrepareSignatureDB(const vector<cv::Point_<T> >& curve, vector<vector<double> >& DB, vector<cv::Point>& DB_params) {
		vector<cv::Point2d> curved; ConvertCurve(curve, curved);
		PrepareSignatureDB(curved, DB, DB_params);
	}

	// prepare curve segment database for generic purposes
	// in particular for comparison between different models
	void PrepareSignatureDB(const vector<cv::Point2d>& curve, vector<vector<double> >& DB, vector<cv::Point2d>& DB_params);

	// prepare curve segment database for comparisons among the shadwos casted by the same model
	// In this case, we would like the min_curve at least half of the given curve.
	void PrepareSignatureDB_For_IntraModel_Comparison
	(const vector<cv::Point2d>& curve, vector<vector<double> >& DB, vector<cv::Point2d>& DB_params);


	//given a match compute transform
	void computeTransform(float a_len, float a_off, float b_len, float b_off, cv::Mat& transform);
	void computeTransform(vector<cv::Point2d>& a_subset, vector<cv::Point2d>& b_subset, cv::Mat& transform);

	//visualize the match
	void visualizeMatching(float a_len, float a_off, float b_len, float b_off);
	//render matching into an image with name "filename"
	void renderMatching(const string& filename, float a_len, float a_off, float b_len, float b_off);
	void drawMatching(cv::Mat& outout, vector<cv::Point2d> & a_subset, vector<cv::Point2d> & b_subset);

	//compute RMSE of the match
	double ComputeRMSE(float a_len, float a_off, float b_len, float b_off);
	double ComputeRMSE(vector<cv::Point2d>& a_subset, vector<cv::Point2d>& b_subset);

	//compute RMSE of the whole curve using the match
	double ComputeWholeRMSE(float a_len, float a_off, float b_len, float b_off);

	//
	//
	// Access functions!
	//
	//

	//get the length of a given curve
	float getLength(const vector<cv::Point2d >& curve);

	void setCurveResampleSize(float size)
	{
		if (size <= 0){
			cerr << "! Warning: CurveArcLengthMatcher::setCurveResampleSize has size < 0" << endl; return;
		}
		RESAMPLE_SIZE = size;
	}
	float getCurveResampleSize() const { return RESAMPLE_SIZE; }

	void setSmallestCurveSegmentSize(float size)
	{
		if (size <= 0){
			cerr << "! Warning: CurveArcLengthMatcher::setSmallestCurveSegmentSize has size < 0" << endl; return;
		}
		SMALLEST_CURVE_SIZE = size;      //smalles curve in the segment DB
	}
	float getSmallestCurveSegmentSize() const { return SMALLEST_CURVE_SIZE; }

	void setLongestCurveSegmentSize(float size)
	{
		if (size <= 0){
			cerr << "! Warning: CurveArcLengthMatcher::setLongestCurveSegmentSize has size < 0" << endl; return;
		}
		LONGEST_CURVE_SIZE = size;
	}
	float getLongestCurveSegmentSize() const { return LONGEST_CURVE_SIZE; }


	void setOffsetStepSize(float size)
	{
		if (size <= 0){
			cerr << "! Warning: CurveArcLengthMatcher::setOffsetStepSize has size < 0" << endl; return;
		}
		OFFSET_STEP = size;               //curve segment DB offset incremental size
	}
	float getOffsetStepSize() const { return OFFSET_STEP; }

	void setCurvatureFilterSize(int size)
	{
		if (size <= 0) return;
		CURVATURE_FILTER_SIZE = size;   //number of matches that will pass through curvature filter
	}
	int getCurvatureFilterSize() const { return CURVATURE_FILTER_SIZE; }

	//load curve from image
	template <typename T>
	void GetCurveForImage(const cv::Mat& filename, vector<cv::Point_<T> >& curve, bool onlyUpper = true, bool getLower = false)
	{
		vector<cv::Point> curve_2i;
		GetCurveForImage(filename, curve_2i, onlyUpper, getLower);
		ConvertCurve(curve_2i, curve);
	}

	//get a subset of curve
	template<typename PT>
	void getSubset(const vector<PT>& curve, float arclen, float offset, vector<PT>& subcurve);

protected:

	bool fileExists(const std::string& filename);



	template<typename T, typename V>
	cv::Mat_<T> ConvertToMat(const vector<vector<V> >& mnt_DB) {
		cv::Mat_<T> mnt_DB_m(mnt_DB.size(), mnt_DB[0].size());
		for (int i = 0; i < mnt_DB.size(); i++) {
			for (int j = 0; j < mnt_DB[i].size(); j++) {
				mnt_DB_m(i, j) = (T)(mnt_DB[i][j]);
			}
		}
		return mnt_DB_m;
	}

	template<typename V>
	cv::Mat_<double> Find2DRigidTransform(const vector<cv::Point_<V> >& a, const vector<cv::Point_<V> >& b,
		cv::Point_<V>* diff = 0, V* angle = 0, V* scale = 0)
	{
		//use PCA to find relational scale
		cv::Mat_<V> P; cv::Mat(a).reshape(1, a.size()).copyTo(P);
		cv::Mat_<V> Q; cv::Mat(b).reshape(1, b.size()).copyTo(Q);
		PCA a_pca(P, cv::Mat(), CV_PCA_DATA_AS_ROW), b_pca(Q, cv::Mat(), CV_PCA_DATA_AS_ROW);

		//disable scaling
		double s = 1; // sqrt(b_pca.eigenvalues.at<V>(0)) / sqrt(a_pca.eigenvalues.at<V>(0));

		//convert to matrices and subtract mean
		//	cv::Mat_<double> P(a.size(),2),Q(b.size(),2);
		Scalar a_m = Scalar(a_pca.mean.at<V>(0), a_pca.mean.at<V>(1));
		Scalar b_m = Scalar(b_pca.mean.at<V>(0), b_pca.mean.at<V>(1));
		//	for (int i=0; i<a.size(); i++) { P(i,0) = a[i].x - a_m[0]; P(i,1) = a[i].y - a_m[1]; }
		//	for (int i=0; i<b.size(); i++) { Q(i,0) = b[i].x - b_m[0]; Q(i,1) = b[i].y - b_m[1]; }
		P -= repeat((cv::Mat_<V>(1, 2) << a_m[0], a_m[1]), P.rows, 1);
		Q -= repeat((cv::Mat_<V>(1, 2) << b_m[0], b_m[1]), Q.rows, 1);

		//    cout << "new mean for a " << mean(P) << "\n";

		////from http://en.wikipedia.org/wiki/Kabsch_algorithm
		//cv::Mat_<double> A = P.t() * Q;
		//SVD svd(A);
		//cv::Mat_<double> C = svd.vt.t() * svd.u.t();
		//double d = (determinant(C) > 0) ? 1 : -1;
		cv::Mat_<double> R = cv::Mat_<double>::eye(cv::Size(2,2));

		cv::Mat_<double> T = (cv::Mat_<double>(3, 3) << 1, 0, b_m[0] / s, 0, 1, b_m[1] / s, 0, 0, 1) *
			(cv::Mat_<double>(3, 3) << s, 0, 0, 0, s, 0, 0, 0, s) *
			(cv::Mat_<double>(3, 3) << R(0, 0), R(0, 1), 0, R(1, 0), R(1, 1), 0, 0, 0, 1) *
			(cv::Mat_<double>(3, 3) << 1, 0, -a_m[0], 0, 1, -a_m[1], 0, 0, 1)
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


private:

	vector<cv::Point2d> source;
	vector<cv::Point2d> target;

	float source_arclength, target_arclength;

	//parameters used in this matcher
	float RESAMPLE_SIZE;            //curve resample size
	float SMALLEST_CURVE_SIZE;      //smallest curve in the segment DB
	float LONGEST_CURVE_SIZE;       //longest curve in the segment DB
	float OFFSET_STEP;              //curve segment DB offset incremental size
	int CURVATURE_FILTER_SIZE;      //number of matches that will pass through curvature filter
};


#endif // ARC_LENGTH_CURVE_SIGNATURE_H_
