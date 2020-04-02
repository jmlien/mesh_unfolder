/*
 *  CurveSignature.cpp
 *  CurveMatching
 *
 *  Created by Roy Shilkrot on 12/7/12.
 *  Copyright 2012 MIT. All rights reserved.
 *
 */

#ifdef _WIN32
#pragma warning( disable : 4018 4244)
#endif

#include "CurveMatching/opencvstd.h"
#include <opencv2/features2d/features2d.hpp>
#include <opencv2/flann/flann.hpp>
#include <opencv2/flann/dist.h>
using namespace cv;

#include "CurveMatching/CurveCSS.h"
#include "CurveMatching/CurveSignature.h"

#include <functional>


//---
// set source of the match
//---
void CurveMatcher::setSource(const string& filename)
{
	vector<Point> a;
	Mat src = imread(filename.c_str());
	if (src.empty()) {
		cerr << "! Error: can't read file " << filename << endl;
		exit(1);
	}

	GetCurveForImage(src, a, false);
	ResampleCurve(a, source, RESAMPLE_SIZE, false);
}

void CurveMatcher::setSource(const vector<Point>& a)
{
	ResampleCurve(a, source, RESAMPLE_SIZE, false);
}

void CurveMatcher::setSource(const vector<Point2d>& a)
{
	if (a.size() != RESAMPLE_SIZE)
		ResampleCurve(a, source, RESAMPLE_SIZE, false);
	else
		source = a;
}

void CurveMatcher::setTarget(const string& filename)
{
	vector<Point> b;
	Mat src2 = imread(filename.c_str());
	if (src2.empty()) {
		cerr << "! Error: can't read file " << filename << endl;
		exit(1);
	}

	GetCurveForImage(src2, b, false);
	ResampleCurve(b, target, RESAMPLE_SIZE, false);
}

void CurveMatcher::setTarget(const vector<Point>& b)
{
	ResampleCurve(b, target, RESAMPLE_SIZE, false);
}

void CurveMatcher::setTarget(const vector<Point2d>& b)
{
	if (b.size() != RESAMPLE_SIZE)
		ResampleCurve(b, target, RESAMPLE_SIZE, false);
	else
		target = b;
}

// prepare curve segment database for comparisons among the shadwos casted by the same model
// In this case, we would like the min_curve at least half of the given curve.
void CurveMatcher::PrepareSignatureDB_For_IntraModel_Comparison
(const vector<cv::Point2d>& curve, vector<vector<double> >& DB, vector<cv::Point>& DB_params)
{

	//automatically determine the value of SMALLEST_CURVE_SIZE
	float SMALLEST_CURVE_SIZE_BK = SMALLEST_CURVE_SIZE;

	SMALLEST_CURVE_SIZE = ((int)floor((RESAMPLE_SIZE / 2) / OFFSET_STEP)) * OFFSET_STEP;
	if (SMALLEST_CURVE_SIZE == 0) SMALLEST_CURVE_SIZE = OFFSET_STEP;

	LONGEST_CURVE_SIZE = RESAMPLE_SIZE;
	//
	PrepareSignatureDB(curve, DB, DB_params);
}

inline void getSubInterval(const vector<Point2d> & curve, int off, int len, vector<Point2d> & sub )
{
	int curve_size = curve.size();
	if(off + len<curve_size)
	{
		sub = vector<Point2d>(curve.begin() + off, curve.begin() + off + len);
	}
	else
	{
		sub = vector<Point2d>(curve.begin() + off, curve.end());
		int remain = off+len-curve_size+1;
		sub.insert(sub.end(),curve.begin(), curve.begin()+remain);
	}
}


// prepare curve segment database for generic purposes
// in particular for comparison between different models
void CurveMatcher::PrepareSignatureDB(const vector<Point2d>& curve_, vector<vector<double> >& DB, vector<Point>& DB_params)
{
	vector<Point2d> curve;
	if (curve_.size() != RESAMPLE_SIZE)
	{
		ResampleCurve(curve_, curve, RESAMPLE_SIZE, true);
	}
	else
	{
		curve = curve_;
	}


	//cout << "curve_ size=" << curve_.size() << endl;

	vector<double> kappa;
	vector<Point2d> small;

	DB.clear(); DB_params.clear();
	int curve_size = curve.size();
	int max_curve_len = min(LONGEST_CURVE_SIZE, curve_size); //curve length cannot exceed curve_size

	for (int len = SMALLEST_CURVE_SIZE; len < max_curve_len - 2; len += 10)
	{
		//int len = RESAMPLE_SIZE;

		//iterate different curve sizes, starting at 20 points
//		cout << "len " << len <<  endl;

		//for (int off = (curve_size - len); off >= 0; )
		for (int off = curve_size-OFFSET_STEP; off >= 0; )
		{
			//int off = 0;

			//iterate segments on A curve
			//this actually has only length (len-1)
			vector<Point2d> small_smooth_input;
			getSubInterval(curve, off, len, small_smooth_input);

			//resample to N points
			ResampleCurve(small_smooth_input, small, RESAMPLE_SIZE, true);

			//compute curvature
			vector<Point2d> small_smooth;
			ComputeCurveCSS(small, kappa, small_smooth, 1.0, true); //smooth and compute the curvature...
			//vector<double> kappa_(kappa.begin()+1,kappa.end()-1);
			vector<double> kappa_(kappa.begin() + 1, kappa.end() - 1);

#if 0
			//cout << "---------------------------------------------------" << endl;
			//for (int i = 0; i < kappa_.size(); i++) cout << kappa_[i] << "\t";
			//cout << endl;

			//if ((off == 79 && len == 110) || (off == 44 && len == 110) )
			{
				//show smooth curve...
				{
					Mat tmptmp(680, 680, CV_8UC3, Scalar::all(0));
					//drawOpenCurve(tmptmp, vector<Point2d>(small_smooth.begin() + 1, small_smooth.end() - 1), kappa_, 2);
					drawOpenCurve(tmptmp, small_smooth_input, kappa, 2);
					//applyColorMap(tmptmp, tmptmp, COLORMAP_HOT);
					imshow("smooth", tmptmp);
					waitKey(0);
				}
			}
			//
#endif

			DB.push_back(kappa_);
			DB_params.push_back(Point(len,off));

			//
			int new_off = off - OFFSET_STEP;
			if (new_off <= 0 && off > 0) off = 0;
			else off = new_off;

		}//end off
	}

	//cout << "! DB size " << DB.size() << endl;
}

inline double L2dist(const vector<double> & a_DB, const vector<double> & b_DB)
{
	double mydist = 0;
	//for (int k = 0; k < a_DB[queryid].size(); k++)
	for (int k = 0; k < a_DB.size(); k++)
	{
		double d = a_DB[k] - b_DB[k];
		mydist += (d*d);
	}

	return sqrt(mydist);
}

//find K closest match between A and b
void CurveMatcher::CompareCurvesUsingSignatureDB
(const vector<cv::Point>& a_DB_params,
const cv::Point& b_DB_params,
const vector<vector<double> >& a_DB,
const vector<double>& b_DB,
cv::DMatch& best_match
)
{
	vector<DMatch> matches;
	vector< pair<double, int> > scores;

	for (int i = 0; i < a_DB.size(); i++)
	{
		double dist=L2dist(a_DB[i],b_DB);
		scores.push_back(make_pair(dist,i));
	}//end for i

	sort(scores.begin(),scores.end());

	//compute RMSD of the smallest...
	int final_size = min(CURVATURE_FILTER_SIZE, (int)scores.size());
	double min_rmse=DBL_MAX;
	int best_match_id=-1;

	for (int i = 0; i < final_size; i++)
	{
		int a_id=scores[i].second;

		int _a_len = a_DB_params[a_id].x;
		int _a_off = a_DB_params[a_id].y;
		int _b_len = b_DB_params.x;
		int _b_off = b_DB_params.y;

		double rmse = ComputeRMSE(_a_len, _a_off, _b_len, _b_off);
		if (rmse < min_rmse)
		{
			best_match_id = a_id;
			min_rmse = rmse;
		}

	}

	//done
	best_match.distance = min_rmse;
	best_match.queryIdx = best_match_id;
	best_match.trainIdx = 0;
}






//find K closest match between A and b
void CurveMatcher::CompareCurvesUsingSignatureDB_curvature_only
(const vector<cv::Point>& a_DB_params,
const cv::Point& b_DB_params,
const vector<vector<double> >& a_DB,
const vector<double>& b_DB,
cv::DMatch& best_match)
{
	double min_score = FLT_MAX;
	int min_id = -1;

	for (int i = 0; i < a_DB.size(); i++)
	{
		double dist = L2dist(a_DB[i], b_DB);
		if (dist < min_score)
		{
			min_score = dist;
			min_id = i;
		}
	}//end for i


	//done
	best_match.distance = min_score;
	best_match.queryIdx = min_id;
	best_match.trainIdx = 0;
}

void CurveMatcher::CompareCurvesUsingSignatureDB
                                   (const vector<Point>& a_DB_params,
								   const vector<Point>& b_DB_params,
								   const vector<vector<double> >& a_DB,
								   const vector<vector<double> >& b_DB,
								   vector<pair<double,DMatch> >& scores_to_matches
								   )
{


	vector<DMatch> matches;


	BFMatcher matcher(NORM_L2);
	Mat_<float> mnt_DB_m = ConvertToMat<float, double>(a_DB); //curvatures of different segments
	Mat_<float> obj_DB_m = ConvertToMat<float, double>(b_DB); //curvatures of different segments
	Mat_<float> dummy;
	matcher.match(mnt_DB_m, obj_DB_m, matches, dummy);

	/* //for some reason, this code is slower than brute force...
	{
		FlannBasedMatcher matcher;
		Mat_<float> mnt_DB_m = ConvertToMat<float,double>(a_DB); //curvatures of different segments
		Mat_<float> obj_DB_m = ConvertToMat<float,double>(b_DB); //curvatures of different segments
		vector<Mat> obj_DB_mv(1,obj_DB_m);
		matcher.add(obj_DB_mv);
		CV_PROFILE(matcher.train();)
		CV_PROFILE(matcher.match( mnt_DB_m, matches );)
	}
	*/

	//DMatch min_match;

	vector<pair<double,int> > scores;
	for( int i = matches.size()-1; i >=0 ; i-- )
	{
		double d = //(max(1.0,200.0 - (double)(a_DB_params[matches[i].queryIdx].x)) +
//					   max(1.0,200.0 - (double)(b_DB_params[matches[i].trainIdx].x))) +
						matches[i].distance;
		if (std::isnan(d))
		{
			cout << "Nan found" << endl;
			continue;
		}
		scores.push_back(make_pair(d, i));
	}
	sort(scores.begin(),scores.end());

	// 110, 14] B = [110, 79]

	/*
	for (int i = 0; i < a_DB.size(); i++)
	{
		if (a_DB_params[i].x != 100 || a_DB_params[i].y!=4) continue;
		for (int j = 0; j < b_DB.size(); j++)
		{
			if (b_DB_params[j].x != 100 || b_DB_params[j].y != 69) continue;
			cout<<"HA-> "<<L2dist(a_DB[i], b_DB[j]) << " => A=" << a_DB_params[i] << " B=" << b_DB_params[j] << endl;
			cout << Mat(a_DB[i]).t() << endl << Mat(b_DB[j]).t() << endl;
		}
	}
	*/


	scores_to_matches.clear();
	int final_size = min(CURVATURE_FILTER_SIZE, (int)scores.size());
	for (int i = 0; i<final_size; i++) { //only look at the top CURVATURE_FILTER_SIZE matches

		int queryid = matches[scores[i].second].queryIdx;
		int trainIdx = matches[scores[i].second].trainIdx;

		//cout << "mydist=" << L2dist(a_DB[queryid], b_DB[trainIdx]) << endl;

		//cout << scores[i].first << " => A=" << a_DB_params[queryid] << " B=" << b_DB_params[trainIdx] << endl;

		//for (int j = 0; j <= 5; j++)
		{
			//int j = 0;
			const DMatch & new_match = matches[scores[i].second];
			//new_match.queryIdx += j;
			//if (new_match.queryIdx < 0 || new_match.queryIdx >= a_DB.size()) continue;
			scores_to_matches.push_back(make_pair(scores[i].first, new_match));
		}
	}

	//min_match = matches[scores.front().second];


	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//double min_dist = scores.front().first;

	////printf("-- Max dist : %f \n", max_dist );
	//printf("-- Min dist : %f, %d(%d,%d) -> %d(%d,%d) \n",
	//	   min_dist,
	//	   min_match.queryIdx,
	//	   a_DB_params[min_match.queryIdx].x,
	//	   a_DB_params[min_match.queryIdx].y,
	//	   min_match.trainIdx,
	//	   b_DB_params[min_match.trainIdx].x,
	//	   b_DB_params[min_match.trainIdx].y
	//	   );
	//
	//cout << Mat(a_DB[min_match.queryIdx]).t() << endl << Mat(b_DB[min_match.trainIdx]).t() << endl;

	/*
	{
		vector<double> a_sig,b_sig;
		a_sig = a_DB[matches[scores.front().second].queryIdx],
		b_sig = b_DB[matches[scores.front().second].trainIdx];
		ShowMathMLCurves(a_sig, b_sig, "curvatures0");

		a_sig = a_DB[matches[scores[1].second].queryIdx];
		b_sig = b_DB[matches[scores[1].second].trainIdx];
		ShowMathMLCurves(a_sig, b_sig, "curvatures1");

		a_sig = a_DB[matches[scores[2].second].queryIdx];
		b_sig = b_DB[matches[scores[2].second].trainIdx];
		ShowMathMLCurves(a_sig, b_sig, "curvatures2");

		a_sig = a_DB[matches[scores[3].second].queryIdx];
		b_sig = b_DB[matches[scores[3].second].trainIdx];
		ShowMathMLCurves(a_sig, b_sig, "curvatures3");

	}
	*/
}



double CurveMatcher::ComputeRMSE(int a_len, int a_off, int b_len, int b_off)
{
	vector<Point2d> a_subset; //(source.begin() + a_off, source.begin() + a_off + a_len);
	vector<Point2d> b_subset; //(target.begin() + b_off, target.begin() + b_off + b_len);
	getSubInterval(source, a_off, a_len, a_subset);
	getSubInterval(target, b_off, b_len, b_subset);

	ResampleCurve(a_subset, a_subset, RESAMPLE_SIZE, true);
	ResampleCurve(b_subset, b_subset, RESAMPLE_SIZE, true);
	return ComputeRMSE(a_subset, b_subset);

}

double CurveMatcher::ComputeRMSE(vector<Point2d>& a_subset, vector<Point2d>& b_subset)
{
	vector<Point2d> ad, bd;
	ConvertCurve(a_subset, ad);
	ConvertCurve(b_subset, bd);
	Mat trans = Find2DRigidTransform(ad, bd);

	vector<Point2d> a_trans;
	cv::transform(ad, a_trans, trans);

	double rmse = 0;
	for (int pt = 0; pt<a_trans.size(); pt++) {
		rmse += norm(a_trans[pt] - bd[pt]);
	}

	return rmse;
}

//compute RMSE of the whole curve using the match
double CurveMatcher::ComputeWholeRMSE(int a_len, int a_off, int b_len, int b_off)
{
	//compute transform and compute RMSE for the aligned segment
	Mat trans;
	vector<Point2d> source_d, target_d;
	ConvertCurve(source, source_d);
	ConvertCurve(target, target_d);

	double rmse = 0;
	{
		vector<Point2d> a_subset; //(source_d.begin() + a_off, source_d.begin() + a_off + a_len);
		vector<Point2d> b_subset; //(target_d.begin() + b_off, target_d.begin() + b_off + b_len);
		getSubInterval(source_d, a_off, a_len, a_subset);
		getSubInterval(target_d, b_off, b_len, b_subset);

		ResampleCurve(a_subset, a_subset, RESAMPLE_SIZE);
		ResampleCurve(b_subset, b_subset, RESAMPLE_SIZE);
		trans = Find2DRigidTransform(a_subset, b_subset);

		vector<Point2d> a_trans;
		cv::transform(a_subset, a_trans, trans);

		for (int pt = 0; pt<a_trans.size(); pt++) {
			rmse += norm(a_trans[pt] - b_subset[pt]);
		}
	}

	//compute RMSE for the rest of the segment

	vector<Point2d> a_rest;// (source_d.begin() + a_off + a_len, source_d.end());
	vector<Point2d> b_rest;// (target_d.begin() + b_off + b_len, target_d.end());
	int a_rest_off = (a_off + a_len) % source_d.size(); int a_rest_len = source_d.size() - a_len;
	int b_rest_off = (b_off + b_len) % target_d.size(); int b_rest_len = target_d.size() - b_len;
	getSubInterval(source_d, a_rest_off, a_rest_len, a_rest);
	getSubInterval(target_d, b_rest_off, b_rest_len, b_rest);
	//a_rest.insert(a_rest.end(), source_d.begin(), source_d.begin() + a_off);
	//b_rest.insert(b_rest.end(), target_d.begin(), target_d.begin() + b_off);

	vector<Point2d> source_p2d, target_p2d;
	ResampleCurve(a_rest, source_p2d, RESAMPLE_SIZE, true);
	ResampleCurve(b_rest, target_p2d, RESAMPLE_SIZE, true);

	vector<Point2d> source_trans;
	cv::transform(source_p2d, source_trans, trans);

	for (int pt = 0; pt<source_trans.size(); pt++) {
		rmse += (norm(source_trans[pt] - target_p2d[pt]));
	}

	return rmse / (2 * RESAMPLE_SIZE); //compute an average
}

void CurveMatcher::CompareCurvesUsingSignatureDB(const vector<Point2d>& a,
	const vector<Point2d>& b,
	int& a_len,
	int& a_off,
	int& b_len,
	int& b_off,
	double& score
	)
{
	vector<Point> a_DB_params, b_DB_params;
	vector<vector<double> > a_DB, b_DB;
	PrepareSignatureDB(a, a_DB, a_DB_params);
	PrepareSignatureDB(b, b_DB, b_DB_params);

	CompareCurvesUsingSignatureDB(a_DB_params, a_DB, b_DB_params, b_DB, a_len, a_off, b_len, b_off, score);
}

void CurveMatcher::CompareCurvesUsingSignatureDB(
	vector<Point>& a_DB_params, vector<vector<double> > & a_DB,
	vector<Point>& b_DB_params, vector<vector<double> > & b_DB,
	int& a_len,
	int& a_off,
	int& b_len,
	int& b_off,
	double& score
	)
{
	vector<pair<double,DMatch> > scores_to_matches;
	CompareCurvesUsingSignatureDB(a_DB_params,b_DB_params,a_DB,b_DB,scores_to_matches);


	//std::cout << "scores_to_matches.size()=" << scores_to_matches.size() << std::endl;

	//re-rank results by RMSE measure after recovering rigid transformation
	DMatch best;
	best.distance = FLT_MAX;
	for (int i=0; i<scores_to_matches.size(); i++)
	{
		int _a_len = a_DB_params[scores_to_matches[i].second.queryIdx].x;
		int _a_off = a_DB_params[scores_to_matches[i].second.queryIdx].y;
		int _b_len = b_DB_params[scores_to_matches[i].second.trainIdx].x;
		int _b_off = b_DB_params[scores_to_matches[i].second.trainIdx].y;

		double rmse = ComputeRMSE(_a_len, _a_off, _b_len, _b_off);
		if (rmse < best.distance)
		{
			best=scores_to_matches[i].second;
			best.distance = rmse;
		}
		//cout << "("<<_a_len<<","<<_a_off<<") -> ("<<_b_len<<","<<_b_off<<")   RMSE: " << rmse << endl;
	}

	//sort(scores_to_matches.begin(), scores_to_matches.end());

	//{
	//	//Show curvatures
	//	vector<double> a_sig,b_sig;
	//	a_sig = a_DB[scores_to_matches.front().second.queryIdx],
	//	b_sig = b_DB[scores_to_matches.front().second.trainIdx];
	//	ShowMathGLCurves(a_sig, b_sig, "curvatures0");
	//}

	a_len = a_DB_params[best.queryIdx].x;
	a_off = a_DB_params[best.queryIdx].y;
	b_len = b_DB_params[best.trainIdx].x;
	b_off = b_DB_params[best.trainIdx].y;
	score = best.distance; // scores_to_matches.front().first;

	//cout << "(" << a_len << "," << a_off << ") -> (" << b_len << "," << b_off << ")   RMSE: " << score << endl;
}


//given a match compute transform
void CurveMatcher::computeTransform(int a_len, int a_off, int b_len, int b_off, Mat& transform)
{
	//Get matched subsets of curves
	vector<Point2d> a_subset(source.begin() + a_off, source.begin() + a_off + a_len);
	vector<Point2d> b_subset(target.begin() + b_off, target.begin() + b_off + b_len);

	//Normalize to equal length
	ResampleCurve(a_subset, a_subset, RESAMPLE_SIZE, true);
	ResampleCurve(b_subset, b_subset, RESAMPLE_SIZE, true);

	computeTransform(a_subset, b_subset, transform);
}

void CurveMatcher::computeTransform(vector<Point2d>& a_subset, vector<Point2d>& b_subset, Mat& transform)
{
	//Prepare the curves for finding the transformation
	vector<Point2f> seq_a_32f, seq_b_32f, seq_a_32f_, seq_b_32f_;

	ConvertCurve(a_subset, seq_a_32f_);
	ConvertCurve(b_subset, seq_b_32f_);

	assert(seq_a_32f_.size() == seq_b_32f_.size());

	seq_a_32f.clear(); seq_b_32f.clear();
	for (int i = 0; i<seq_a_32f_.size(); i++) {
		//		if(i%2 == 0) { // you can use only part of the points to find the transformation
		seq_a_32f.push_back(seq_a_32f_[i]);
		seq_b_32f.push_back(seq_b_32f_[i]);
		//		}
	}
	assert(seq_a_32f.size() == seq_b_32f.size()); //just making sure

	vector<Point2d> seq_a_trans(source.size());

	//Find the fitting transformation
	//	Mat affineT = estimateRigidTransform(seq_a_32f,seq_b_32f,false); //may wanna use Affine here..
	transform = Find2DRigidTransform(seq_a_32f, seq_b_32f);
}

void CurveMatcher::visualizeMatching(int a_len, int a_off, int b_len, int b_off)
{
	//Get matched subsets of curves
	vector<Point2d> a_subset; //(source.begin() + a_off, source.begin() + a_off + a_len);
	vector<Point2d> b_subset; //(target.begin() + b_off, target.begin() + b_off + b_len);
	getSubInterval(source, a_off, a_len, a_subset);
	getSubInterval(target, b_off, b_len, b_subset);

	//vector<Point2d> a_subset(source.begin() + a_off, source.begin() + a_off + a_len);
	//vector<Point2d> b_subset(target.begin() + b_off, target.begin() + b_off + b_len);

	//Normalize to equal length
	ResampleCurve(a_subset, a_subset, RESAMPLE_SIZE, true);
	ResampleCurve(b_subset, b_subset, RESAMPLE_SIZE, true);

	//
	Mat outout(600, 600, CV_8UC3, Scalar::all(0));
	drawMatching(outout, a_subset, b_subset);

//#ifndef __APPLE__
	imshow("outout", outout);
	waitKey();
//#endif
}

void CurveMatcher::renderMatching(const string& filename, int a_len, int a_off, int b_len, int b_off)
{
	//Get matched subsets of curves
	vector<Point2d> a_subset; //(source.begin() + a_off, source.begin() + a_off + a_len);
	vector<Point2d> b_subset; //(target.begin() + b_off, target.begin() + b_off + b_len);
	getSubInterval(source, a_off, a_len, a_subset);
	getSubInterval(target, b_off, b_len, b_subset);
	//vector<Point2d> a_subset(source.begin() + a_off, source.begin() + a_off + a_len);
	//vector<Point2d> b_subset(target.begin() + b_off, target.begin() + b_off + b_len);

	//Normalize to equal length
	ResampleCurve(a_subset, a_subset, RESAMPLE_SIZE, true);
	ResampleCurve(b_subset, b_subset, RESAMPLE_SIZE, true);

	//
	Mat outout(600, 600, CV_8UC3, Scalar::all(0));
	drawMatching(outout, a_subset, b_subset);
	imwrite(filename.c_str(), outout);
}

void CurveMatcher::drawMatching(Mat& outout, vector<Point2d> & a_subset, vector<Point2d> & b_subset)
{
	vector<Point2d> a_p2d;
	vector<Point2d> b_p2d;
	ConvertCurve(source, a_p2d); //convert to point
	ConvertCurve(target, b_p2d); //convert to point

	//Visualize the original and target

	{
		// //draw small original
		// vector<Point2d> tmp_curve;
		// cv::transform(a_p2d, tmp_curve, getRotationMatrix2D(Point2f(0, 0), 0, .5));
		// Mat tmp_curve_m(tmp_curve); tmp_curve_m += Scalar(5, 0);
		// drawOpenCurve(outout, tmp_curve, Scalar(255), 1);
		//
		// //draw small matched subset of original
		// ConvertCurve(a_subset, tmp_curve);
		// cv::transform(tmp_curve, tmp_curve, getRotationMatrix2D(Point2f(0, 0), 0, .5));
		// Mat tmp_curve_m1(tmp_curve); tmp_curve_m1 += Scalar(5, 0);
		// drawOpenCurve(outout, tmp_curve, Scalar(255, 255), 2);
		//
		// //draw small target
		// ConvertCurve(target, tmp_curve);
		// cv::transform(tmp_curve, tmp_curve, getRotationMatrix2D(Point2f(0, 0), 0, .5));
		// Mat tmp_curve_m2(tmp_curve); tmp_curve_m2 += Scalar(outout.cols - 300, 0);
		// drawOpenCurve(outout, tmp_curve, Scalar(255, 0, 255), 1);
		//
		// //draw small matched subset of target
		// ConvertCurve(b_subset, tmp_curve);
		// cv::transform(tmp_curve, tmp_curve, getRotationMatrix2D(Point2f(0, 0), 0, .5));
		// Mat tmp_curve_m3(tmp_curve); tmp_curve_m3 += Scalar(outout.cols - 300, 0);
		// drawOpenCurve(outout, tmp_curve, Scalar(255, 175, 255), 2);

		//draw big target
		drawOpenCurve(outout, target, Scalar(0, 0, 255), 1);
		//draw big matched subset of target
		drawOpenCurve(outout, b_subset, Scalar(0, 255, 255), 2);
	}


	//Prepare the curves for finding the transformation
	Mat trans;
	computeTransform(a_subset, b_subset, trans);
	//cout << trans;

	vector<Point2d> seq_a_trans(source.size());
	cv::transform(a_p2d, seq_a_trans, trans);

	vector<Point2d> a_subset_trans(source.size());
	vector<Point2d> tmp_curve;
	ConvertCurve(a_subset, tmp_curve);
	cv::transform(tmp_curve, a_subset_trans, trans);

	//draw the result matching : the complete original curve as matched to the target
	drawOpenCurve(outout, seq_a_trans,    Scalar(0, 255, 0), 1);
	drawOpenCurve(outout, a_subset_trans, Scalar(125, 255, 125), 2);

	//May want to visualize point-by-point matching
	//	cv::transform(seq_a_32f,seq_a_32f,trans);
	//	for (int i=0; i<seq_a_32f.size(); i++) {
	//		line(outout, seq_a_32f[i], seq_b_32f[i], Scalar(0,0,255), 1);
	//	}
}


Mat_<double> CurveMatcher::GetSmithWatermanHMatrix(const vector<pair<char, int> >& a, const vector<pair<char, int> >& b) {
	int M = a.size();
	int N = b.size();

	//Smith-Waterman
	Mat_<double> H(M + 1, N + 1, 0.0);
	for (int i = 1; i <= M; i++) {
		for (int j = 1; j <= N; j++) {
			vector<double> v(4, 0.0);
			v[1] = H(i - 1, j - 1) + ((a[i - 1].first == b[j - 1].first) ? 2.0 : -1.0);
			v[2] = H(i - 1, j) - 1.0;
			v[3] = H(i, j - 1) - 1.0;
			H(i, j) = *(max_element(v.begin(), v.end()));
		}
	}
	//	cout << H << endl;
	return H;
}

//JML: this is interesting keep it
/* original Smith Waterman algorithm */
double CurveMatcher::MatchSmithWaterman(const vector<pair<char, int> >& a, const vector<pair<char, int> >& b, vector<Point>& matching)
{
	vector<Point> traceback;
	Mat_<double> H = GetSmithWatermanHMatrix(a, b);
	Point maxp; double maxval;
	minMaxLoc(H, NULL, &maxval, NULL, &maxp);
	vector<char> step;
	while (H(maxp.y, maxp.x) != 0) {
		//				cout << "H(maxp.y-1,maxp.x-1) > H(maxp.y,maxp.x-1)" << H(maxp.y-1,maxp.x-1) << " > " << H(maxp.y,maxp.x-1) << endl;
		if (H(maxp.y - 1, maxp.x - 1) > H(maxp.y, maxp.x - 1) &&
			H(maxp.y - 1, maxp.x - 1) > H(maxp.y - 1, maxp.x))
		{
			traceback.push_back(maxp);
			maxp = maxp - Point(1, 1);
			step.push_back('a');
		}
		else
			if (H(maxp.y - 1, maxp.x) > H(maxp.y - 1, maxp.x - 1) &&
				H(maxp.y - 1, maxp.x) > H(maxp.y, maxp.x - 1))
			{
				traceback.push_back(maxp);
				maxp.y--;
				step.push_back('d');
			}
			else
				if (H(maxp.y, maxp.x - 1) > H(maxp.y - 1, maxp.x - 1) &&
					H(maxp.y, maxp.x - 1) > H(maxp.y - 1, maxp.x))
				{
					traceback.push_back(maxp);
					maxp.x--;
					step.push_back('i');
				}
				else {
					//default - go back on both
					traceback.push_back(maxp);
					maxp = maxp - Point(1, 1);
					step.push_back('a');
				}
	}
	for (vector<Point>::reverse_iterator it = traceback.rbegin();
		it != traceback.rend() - 1;
		++it)
	{
		if ((*it).y != (*(it + 1)).y && (*it).x != (*(it + 1)).x)
			matching.push_back(Point((*it).y, (*it).x));
	}
	for (vector<Point>::reverse_iterator it = traceback.rbegin();
		it != traceback.rend();
		++it)
	{
		if (it == traceback.rend())
			cout << a[(*it).y].first;
		else {
			if ((*it).y == (*(it + 1)).y)
				cout << "-";
			else {
				cout << a[(*it).y].first;
			}
		}
	}
	cout << endl;
	for (vector<Point>::reverse_iterator it = traceback.rbegin();
		it != traceback.rend();
		++it)
	{
		if (it == traceback.rend())
			cout << b[(*it).x].first;
		else {
			if ((*it).x == (*(it + 1)).x)
				cout << "-";
			else
				cout << b[(*it).x].first;
		}
	}
	cout << endl;
	for (int k = 0; k<step.size(); k++) {
		cout << step[k];
	}
	cout << endl;
	return maxval;
}


bool CurveMatcher::fileExists(const std::string& filename)
{
	struct stat buf;
	if (stat(filename.c_str(), &buf) != -1)
	{
		return true;
	}
	return false;
}

void CurveMatcher::GetCurveForImage(const Mat& filename, vector<Point>& whole, vector<Point>& curve_upper, vector<Point>& curve_lower)
{
	assert(!filename.empty());
	Mat tmp; filename.copyTo(tmp);
	Mat gray;
	if (tmp.type() == CV_8UC3)
		cvtColor(tmp, gray, COLOR_BGR2GRAY);
	else if (tmp.type() == CV_8UC1)
		gray = tmp;
	else
		CV_Error(-1, "GetCurveForImage: unsupported image format"); //, __FILE__, __LINE__);


	threshold(gray, gray, 128, 255, THRESH_BINARY);
	//	imshow("input",gray);

	//negate the color
	bitwise_not(gray, gray);

	vector<vector<Point> > contours;
	findContours(gray, contours, RETR_EXTERNAL, CHAIN_APPROX_SIMPLE);
	if (contours.size() <= 0) return;

	vector<Point> upperCurve = contours[0];
	if (upperCurve.size() <= 50) {
		return;
	}

	//find minimal and maximal X coord
	vector<double> x, y;
	PolyLineSplit(contours[0], x, y);
	Point minxp, maxxp;
	minMaxLoc(x, 0, 0, &minxp, &maxxp);
	int minx = minxp.x, maxx = maxxp.x;
	if (minx > maxx) swap(minx, maxx);

	//take lower and upper halves of the curve
	vector<Point> upper, lower;
	upper.insert(upper.begin(), contours[0].begin() + minx, contours[0].begin() + maxx);
	lower.insert(lower.begin(), contours[0].begin() + maxx, contours[0].end());
	lower.insert(lower.end(), contours[0].begin(), contours[0].begin() + minx);

	//test which is really the upper part, by looking at the y-coord of the mid point

	if (lower[lower.size() / 2].y <= upper[upper.size() / 2].y) {
		curve_upper = lower;
		curve_lower = upper;
	}
	else {
		curve_upper = upper;
		curve_lower = lower;
	}

	//make sure it goes left-to-right
	if (curve_upper.front().x > curve_upper.back().x) { //hmmm, need to flip
		reverse(curve_upper.begin(), curve_upper.end());
	}

	whole.clear();
	whole.insert(whole.begin(), curve_upper.rbegin(), curve_upper.rend());
	whole.insert(whole.begin(), curve_lower.begin(), curve_lower.end());
}

void CurveMatcher::GetCurveForImage(const Mat& filename, vector<Point>& curve, bool onlyUpper, bool getLower)
{
	vector<Point> whole, upper, lower;
	GetCurveForImage(filename, whole, upper, lower);
	if (onlyUpper) {
		if (getLower)
			curve = lower;
		else
			curve = upper;
	}
	else {
		curve = whole;
	}
}
