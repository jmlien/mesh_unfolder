
#ifdef _WIN32
#pragma warning( disable : 4018 4244)
#endif

#include "CurveMatching/ArcLengthCurveSignature.h"
#include "CurveMatching/opencvstd.h"
#include <opencv2/features2d/features2d.hpp>
#include <opencv2/flann/flann.hpp>
#include <opencv2/flann/dist.h>
using namespace cv;

#include "CurveMatching/CurveCSS.h"

//-----------------------------------------------------------------------------
//
//
// debug parameters in this file
//
//
//-----------------------------------------------------------------------------

#define SHOW_RMSE_MATCHING_ANIMATION    00

extern bool ARC_LEN_DEBUG;

#include <functional>

float getArcLength(const vector<cv::Point2d>& curve)
{
	if (curve.size()<=1) return 0;
	auto pre = curve.begin();
	float len = 0;
	for (auto i = ++curve.begin(); i != curve.end(); i++)
	{
		double d = norm(*i - *pre);
		len += d;
		pre = i;
	}
	return len;
}

//---
// set source of the match
//---
void CurveArcLengthMatcher::setSource(const string& filename)
{
	vector<Point> a;
	Mat src = imread(filename.c_str());
	if (src.empty()) {
		cerr << "! Error: can't read file " << filename << endl;
		exit(1);
	}

	GetCurveForImage(src, a, false);
	//vector<cv::Point2d> tmp;
	ConvertCurve(a, source);

	source_arclength = getArcLength(source);
	//ResampleCurve(tmp, source, RESAMPLE_SIZE, false);
}

void CurveArcLengthMatcher::setSource(const vector<cv::Point2d>& a)
{
	//if (a.size() != RESAMPLE_SIZE)
	//	ResampleCurve(a, source, RESAMPLE_SIZE, false);
	//else
		source = a;
		source_arclength = getArcLength(source);
}

void CurveArcLengthMatcher::setTarget(const string& filename)
{
	vector<Point> b;
	Mat src2 = imread(filename.c_str());
	if (src2.empty()) {
		cerr << "! Error: can't read file " << filename << endl;
		exit(1);
	}

	GetCurveForImage(src2, b, false);
	//vector<cv::Point2d> tmp;
	ConvertCurve(b, target);
	target_arclength = getArcLength(target);
	//ResampleCurve(tmp, target, RESAMPLE_SIZE, false);
}

void CurveArcLengthMatcher::setTarget(const vector<cv::Point2d>& b)
{
	//if (b.size() != RESAMPLE_SIZE)
	//	ResampleCurve(b, target, RESAMPLE_SIZE, false);
	//else
	target = b;
	target_arclength = getArcLength(target);
}

//
// get a subset of curve based the the offset and the total arc length of the subset
//
template<typename PT>
void CurveArcLengthMatcher::getSubset(const vector<PT>& curve, float arclen, float offset, vector<PT>& subcurve)
{
	//get a subset of curve
	float len = 0;
	auto pre = curve.rbegin();
	float start = offset;
	float end = (arclen>0)?(offset + arclen):FLT_MAX;

	for (auto i = ++curve.rbegin(); i != curve.rend();i++)
	{
		double d=cv::norm(*i - *pre);
		double new_len = len + d;
		if (len <= start && start <= new_len) //find the offset point
		{
			float t = ((start - len) / (d));
			PT pt = t*(*i) + (1 - t)*(*pre);
			subcurve.push_back(pt);
			//subcurve.push_back(*i);
		}

		if (new_len >= start && new_len <= end)
		{
			subcurve.push_back(*i);
		}
		else if (new_len>end) //find the end
		{
			float t = ((end - len) / (d));
			PT pt = t*(*i) + (1 - t)*(*pre);
			subcurve.push_back(pt);
			break;
		}

		len += d;
		pre = i;
	}

	if (subcurve.empty())
	{
		cerr << "! Warning: subcurve is empty len=" << len << " arclen="<<arclen<<"  offset="<<offset<<" end="<<offset+arclen << " arclen(curve)="
			<< getArcLength(curve) << endl;
	}
}

// prepare curve segment database for comparisons among the shadwos casted by the same model
// In this case, we would like the min_curve at least half of the given curve.
void CurveArcLengthMatcher::PrepareSignatureDB_For_IntraModel_Comparison
(const vector<cv::Point2d>& curve_, vector<vector<double> >& DB, vector<cv::Point2d>& DB_params)
{
	vector<cv::Point2d> curve;
	if (curve_.size() != RESAMPLE_SIZE) {
		ResampleCurve(curve_, curve, RESAMPLE_SIZE, true);
	}
	else {
		curve = curve_;
	}

	//automatically determine the value of SMALLEST_CURVE_SIZE
	float curve_arclen = getArcLength(curve);
	float SMALLEST_CURVE_SIZE_BK = SMALLEST_CURVE_SIZE;

	SMALLEST_CURVE_SIZE = ((int)floor((curve_arclen / 2) / OFFSET_STEP)) * OFFSET_STEP;
	if (SMALLEST_CURVE_SIZE == 0) SMALLEST_CURVE_SIZE = OFFSET_STEP;

	//cout << "!!!!!! # curve_arclen=" << curve_arclen << " OFFSET_STEP=" << OFFSET_STEP << " SMALLEST_CURVE_SIZE=" << SMALLEST_CURVE_SIZE<< endl;

	LONGEST_CURVE_SIZE = curve_arclen;
	//
	PrepareSignatureDB(curve_, DB, DB_params);
}

//
// prepare curve segment database for generic purposes
// in particular for comparison between different models
//
void CurveArcLengthMatcher::PrepareSignatureDB(const vector<cv::Point2d>& curve_, vector<vector<double> >& DB, vector<cv::Point2d>& DB_params)
{
	vector<cv::Point2d> curve=curve_;
	vector<double> kappa;
	vector<cv::Point2d> small;

	DB.clear(); DB_params.clear();
	float curve_arclen = getArcLength(curve);

	//cout << "CUVE ARC LENGTH = " << curve_arclen << endl;

	if (curve_arclen < SMALLEST_CURVE_SIZE)
	{
		cout << "SMALLEST_CURVE_SIZE=" << SMALLEST_CURVE_SIZE << " curve_arclen=" << curve_arclen << endl;
	}

	float max_curve_len = min(LONGEST_CURVE_SIZE, curve_arclen); //LONGEST_CURVE_SIZE cannot be larger than curve_arclen

	//JML: TODO: the value 10 below should be a parameter instead of a constant
	for (float len = SMALLEST_CURVE_SIZE; len <= max_curve_len; len += 10) //len is arc length
	{
		//int len = RESAMPLE_SIZE;

		//iterate different curve sizes, starting at 20 points
		//cout << "CURVE SEG LEN " << len <<  endl;

		for (float off = (curve_arclen - len); off >= 0;)
		{
			//int off = 0;

			//JML: TODO: I think we should try getSubset(curve_, ...) below instead of getSubset(curve,...)
			//iterate segments on A curve
			//this actually has only length (len-1)
			vector<cv::Point2d> small_smooth_input;
			getSubset(curve, len, off, small_smooth_input);

			//resample to N points
			ResampleCurve(small_smooth_input, small, RESAMPLE_SIZE, true);

			//compute curvature
			vector<cv::Point2d> small_smooth;
			ComputeCurveCSS(small, kappa, small_smooth, 5.0, true); //smooth and compute the curvature...
			vector<double> kappa_(kappa.begin()+1,kappa.end()-1);

#if 0
			//cout << "---------------------------------------------------" << endl;
			//for (int i = 0; i < kappa_.size(); i++) cout << kappa_[i] << "\t";
			//cout << endl;

			if ((off == 79 && len == 110) || (off == 44 && len == 110) )
			{
				//show smooth curve...
				{
					Mat tmptmp(680, 680, CV_8UC3, Scalar::all(0));
					//drawOpenCurve(tmptmp, vector<cv::Point2d>(small_smooth.begin() + 1, small_smooth.end() - 1), kappa_, 2);
					drawOpenCurve(tmptmp, small_smooth_input, kappa, 2);
					//applyColorMap(tmptmp, tmptmp, COLORMAP_HOT);
					imshow("smooth", tmptmp);
					waitKey(0);
				}
			}
			//
#endif

			DB.push_back(kappa_);
			DB_params.push_back(cv::Point2d(len,off));

			//
			float new_off = off - OFFSET_STEP;
			if (new_off <= 0 && off > 0) off = 0;
			else off = new_off;

		}//end off
	}

	//cout << "! DB size " << DB.size() << endl;
}

float CurveArcLengthMatcher::getLength(const vector<cv::Point2d >& curve)
{
	return getArcLength(curve);
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
void CurveArcLengthMatcher::CompareCurvesUsingSignatureDB
(const vector<cv::Point2d>& a_DB_params,
const cv::Point2d& b_DB_params,
const vector<vector<double> >& a_DB,
const vector<double>& b_DB,
DMatch& best_match
)
{
	//
	//Note: Only compare curve segment with the same arc length
	//
	vector<DMatch> matches;
	vector< pair<double, int> > scores;

	for (int i = 0; i < a_DB.size(); i++)
	{
		if (a_DB_params[i].x != b_DB_params.x) continue; //the arc length is different
		double dist=L2dist(a_DB[i],b_DB);
		scores.push_back(make_pair(dist,i));
	}//end for i

	if (scores.empty())
	{
		//handle the case that no match is found
		best_match.distance = FLT_MAX;
		best_match.queryIdx = -1;
		return;
	}

	//now, there are some matches, find the best one
	sort(scores.begin(),scores.end());

	//compute RMSD of the smallest...
	int final_size = min(CURVATURE_FILTER_SIZE, (int)scores.size());
	double min_rmse=DBL_MAX;
	int best_match_id=0;

	for (int i = 0; i < final_size; i++)
	{
		int a_id=scores[i].second;

		float _a_len = a_DB_params[a_id].x;
		float _a_off = a_DB_params[a_id].y;
		float _b_len = b_DB_params.x;
		float _b_off = b_DB_params.y;

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
void CurveArcLengthMatcher::CompareCurvesUsingSignatureDB_curvature_only
(const vector<cv::Point2d>& a_DB_params,
const cv::Point2d& b_DB_params,
const vector<vector<double> >& a_DB,
const vector<double>& b_DB,
DMatch& best_match
)
{
	vector<DMatch> matches;
	vector< pair<double, int> > scores;

	for (int i = 0; i < a_DB.size(); i++)
	{
		if (a_DB_params[i].x != b_DB_params.x) continue; //the arc length is different
		double dist = L2dist(a_DB[i], b_DB);
		scores.push_back(make_pair(dist, i));
	}//end for i

	if (scores.empty())
	{
		//handle the case that no match is found
		best_match.distance = FLT_MAX;
		best_match.queryIdx = -1;
		return;
	}

	//now, there are some matches, find the best one
	sort(scores.begin(), scores.end());

	//compute RMSD of the smallest...
	int final_size = min(CURVATURE_FILTER_SIZE, (int)scores.size());
	double min_rmse = DBL_MAX;
	int best_match_id = 0;

	for (int i = 0; i < final_size; i++)
	{
		int a_id = scores[i].second;

		float _a_len = a_DB_params[a_id].x;
		float _a_off = a_DB_params[a_id].y;
		float _b_len = b_DB_params.x;
		float _b_off = b_DB_params.y;

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

//
// only compare between curve segments with the same arc length....
//

void CurveArcLengthMatcher::CompareCurvesUsingSignatureDB
(
const vector<cv::Point2d>& a_DB_params,
const vector<cv::Point2d>& b_DB_params,
const vector<vector<double> >& a_DB,
const vector<vector<double> >& b_DB,
vector<pair<double,DMatch> >& scores_to_matches
)
{
	//mapping from curve segment length to DB indice
	map<double, vector<int> > map_len_offset_a;
	map<double, vector<int> > map_len_offset_b;

	{
		int index = 0;
		for (const auto & i : a_DB_params)
			map_len_offset_a[i.x].push_back(index++);

		index = 0;
		for (const auto & i : b_DB_params)
			map_len_offset_b[i.x].push_back(index++);
	}
	//

	//find matches for the curves with the same length
	vector<DMatch> matches;
	for (const auto & i : map_len_offset_a)
	{
		if (map_len_offset_b.find(i.first) == map_len_offset_b.end()) continue; // nothing to match to...
		const vector<int>& a_list = i.second;
		const vector<int>& b_list = map_len_offset_b[i.first];

		//collect the subset and perform matching
		vector<cv::Point2d> a_DB_params_subset;
		vector<cv::Point2d> b_DB_params_subset;
		vector<vector<double> > a_DB_subset;
		vector<vector<double> > b_DB_subset;

		//
		for (auto a : a_list)
		{
			a_DB_params_subset.push_back(a_DB_params[a]);
			a_DB_subset.push_back(a_DB[a]);
		}

		for (auto b : b_list)
		{
			b_DB_params_subset.push_back(b_DB_params[b]);
			b_DB_subset.push_back(b_DB[b]);
		}

		//perform maching here


		//
		vector<DMatch> local_matches;
		BFMatcher matcher(NORM_L2);
		Mat_<float> mnt_DB_m = ConvertToMat<float, double>(a_DB_subset); //curvatures of different segments
		Mat_<float> obj_DB_m = ConvertToMat<float, double>(b_DB_subset); //curvatures of different segments
		Mat_<float> dummy;
		matcher.match(mnt_DB_m, obj_DB_m, local_matches, dummy);

		//insert the local matches to the global list
		for (auto& m : local_matches)
		{
			m.queryIdx = a_list[m.queryIdx];
			m.trainIdx = b_list[m.trainIdx];
			matches.push_back(m);
		}
		//matches.insert(matches.end(), local_matches.begin(), local_matches.end());
	}



	//
	vector<pair<double,int> > scores;
	for( int i = matches.size()-1; i >=0 ; i-- )
	{
		double d = matches[i].distance;
		if (isnan(d))
		{
			cout << "Nan found" << endl;
			continue;
		}
		scores.push_back(make_pair(d, i));
	}
	sort(scores.begin(),scores.end());

	scores_to_matches.clear();
	int final_size = min(CURVATURE_FILTER_SIZE, (int)scores.size());
	for (int i = 0; i<final_size; i++)  //only look at the top CURVATURE_FILTER_SIZE matches
	{
		//int queryid = matches[scores[i].second].queryIdx;
		//int trainIdx = matches[scores[i].second].trainIdx;

		const DMatch & new_match = matches[scores[i].second];
		scores_to_matches.push_back(make_pair(scores[i].first, new_match));
	}

	//min_match = matches[scores.front().second];


	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//double min_dist = scores.front().first;
}



double CurveArcLengthMatcher::ComputeRMSE(float a_len, float a_off, float b_len, float b_off)
{
	vector<cv::Point2d> a_subset, b_subset;

	getSubset(source, a_len, a_off, a_subset);
	getSubset(target, b_len, b_off, b_subset);

	ResampleCurve(a_subset, a_subset, RESAMPLE_SIZE, true);
	ResampleCurve(b_subset, b_subset, RESAMPLE_SIZE, true);
	return ComputeRMSE(a_subset, b_subset);

}

double CurveArcLengthMatcher::ComputeRMSE(vector<cv::Point2d>& a_subset, vector<cv::Point2d>& b_subset)
{
	Mat trans = Find2DRigidTransform(a_subset, b_subset);

	vector<cv::Point2d> a_trans;
	cv::transform(a_subset, a_trans, trans);
	const auto size = a_trans.size();
	double rmse = 0;
	for (int pt = 0; pt<size; pt++) {
		const auto d = a_trans[pt] - b_subset[pt];
		rmse += (d.x*d.x+d.y*d.y);
	}

	return sqrt(rmse / size);
}


//compute RMSE of the whole curve using the match
double CurveArcLengthMatcher::ComputeWholeRMSE(float a_len, float a_off, float b_len, float b_off)
{
	//compute transform and compute RMSE for the aligned segment
	Mat trans;

	{
		vector<cv::Point2d> a_subset, b_subset;
		getSubset(source, a_len, a_off, a_subset);
		getSubset(target, b_len, b_off, b_subset);
		ResampleCurve(a_subset, a_subset, RESAMPLE_SIZE, true);
		ResampleCurve(b_subset, b_subset, RESAMPLE_SIZE, true);
		trans = Find2DRigidTransform(a_subset, b_subset);
	}

	double rmse = 0;

	vector<cv::Point2d> a_subset, b_subset;
	getSubset(source, -1, a_off, a_subset);   //between (a_off+a_len) and the end
	getSubset(target, -1, b_off, b_subset); //between (b_off+b_len) and the end
	if (a_off != 0) getSubset(source, a_off, 0, a_subset);                               //between 0 and a_off
	if (b_off != 0) getSubset(target, b_off, 0, b_subset);


	ResampleCurve(a_subset, a_subset, RESAMPLE_SIZE * 2, false);
	ResampleCurve(b_subset, b_subset, RESAMPLE_SIZE * 2, false);

	//cout << "a_subset 2=";
	//for (int i = 0; i < a_subset.size(); i++)
	//{
	//	cout << "(" << a_subset[i].x << "," << a_subset[i].y << "), ";
	//}
	//cout << endl << endl << endl;

	//cout << "b_subset 2=";
	//for (int i = 0; i < b_subset.size(); i++)
	//{
	//	cout << "(" << b_subset[i].x << "," << b_subset[i].y << "), ";
	//}
	//cout << endl << endl << endl;

	//cout << "source=";
	//for (int i = 0; i < source.size(); i++)
	//{
	//	cout << "(" << source[i].x << "," << source[i].y << "), ";
	//}
	//cout << endl << endl << endl;

	//cout << "target=";
	//for (int i = 0; i < target.size(); i++)
	//{
	//	cout << "(" << target[i].x << "," << target[i].y << "), ";
	//}
	//cout << endl << endl << endl;


	vector<cv::Point2d> a_trans;
	cv::transform(a_subset, a_trans, trans);

	for (int pt = 0; pt<a_trans.size(); pt++) {
		const auto d = a_trans[pt] - b_subset[pt];
		rmse += (d.x*d.x + d.y*d.y);
	}

	//for debugging only
	//show an animation of how points are mapping to each other
#if SHOW_RMSE_MATCHING_ANIMATION
	for (int i = 0; i < RESAMPLE_SIZE * 2; i++)
	{

		cv::Mat outout(1000, 1000, CV_8UC3, cv::Scalar::all(255));
		drawOpenCurve(outout, a_trans, Scalar(255, 0, 0), 5, i);
		drawOpenCurve(outout, b_subset, Scalar(0, 0, 255), 1, i);
		imshow("no", outout);
		cv::waitKey(10);
	}
#endif

	return sqrt(rmse / (2 * RESAMPLE_SIZE)); //compute an average
}

void CurveArcLengthMatcher::CompareCurvesUsingSignatureDB(const vector<cv::Point2d>& a,
	const vector<cv::Point2d>& b,
	float& a_len,
	float& a_off,
	float& b_len,
	float& b_off,
	double& score
	)
{
	vector<cv::Point2d> a_DB_params, b_DB_params;
	vector<vector<double> > a_DB, b_DB;
	PrepareSignatureDB(a, a_DB, a_DB_params);
	PrepareSignatureDB(b, b_DB, b_DB_params);

	CompareCurvesUsingSignatureDB(a_DB_params, a_DB, b_DB_params, b_DB, a_len, a_off, b_len, b_off, score);
}

void CurveArcLengthMatcher::CompareCurvesUsingSignatureDB
(
	vector<cv::Point2d>& a_DB_params, vector<vector<double> > & a_DB,
	vector<cv::Point2d>& b_DB_params, vector<vector<double> > & b_DB,
	float& a_len,
	float& a_off,
	float& b_len,
	float& b_off,
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
		float _a_len = a_DB_params[scores_to_matches[i].second.queryIdx].x;
		float _a_off = a_DB_params[scores_to_matches[i].second.queryIdx].y;
		float _b_len = b_DB_params[scores_to_matches[i].second.trainIdx].x;
		float _b_off = b_DB_params[scores_to_matches[i].second.trainIdx].y;

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

	if (best.distance == FLT_MAX) //no match found...
	{
		score = DBL_MAX;
		return;
	}

	a_len = a_DB_params[best.queryIdx].x;
	a_off = a_DB_params[best.queryIdx].y;
	b_len = b_DB_params[best.trainIdx].x;
	b_off = b_DB_params[best.trainIdx].y;
	score = best.distance; // scores_to_matches.front().first;

	//cout << "(" << a_len << "," << a_off << ") -> (" << b_len << "," << b_off << ")   RMSE: " << score << endl;
}


//given a match compute transform
void CurveArcLengthMatcher::computeTransform(float a_len, float a_off, float b_len, float b_off, Mat& transform)
{
	vector<cv::Point2d> a_subset, b_subset; // (source.begin() + a_off, source.begin() + a_off + a_len);
	//vector<cv::Point2d> b_subset(target.begin() + b_off, target.begin() + b_off + b_len);

	getSubset(source, a_len, a_off, a_subset);
	getSubset(target, b_len, b_off, b_subset);

	ResampleCurve(a_subset, a_subset, RESAMPLE_SIZE, true);
	ResampleCurve(b_subset, b_subset, RESAMPLE_SIZE, true);

	computeTransform(a_subset, b_subset, transform);
}

//
// This can be simplified...
//
void CurveArcLengthMatcher::computeTransform(vector<cv::Point2d>& seq_a_32f_, vector<cv::Point2d>& seq_b_32f_, Mat& transform)
{
	//Prepare the curves for finding the transformation
	vector<Point2f> seq_a_32f, seq_b_32f; // , seq_a_32f_, seq_b_32f_;

	//ConvertCurve(a_subset, seq_a_32f_);
	//ConvertCurve(b_subset, seq_b_32f_);

	assert(seq_a_32f_.size() == seq_b_32f_.size());

	seq_a_32f.clear(); seq_b_32f.clear();
	for (int i = 0; i<seq_a_32f_.size(); i++) {
		//		if(i%2 == 0) { // you can use only part of the points to find the transformation
		seq_a_32f.push_back(seq_a_32f_[i]);
		seq_b_32f.push_back(seq_b_32f_[i]);
		//		}
	}
	assert(seq_a_32f.size() == seq_b_32f.size()); //just making sure

	vector<cv::Point2d> seq_a_trans(source.size());

	//Find the fitting transformation
	//	Mat affineT = estimateRigidTransform(seq_a_32f,seq_b_32f,false); //may wanna use Affine here..
	transform = Find2DRigidTransform(seq_a_32f, seq_b_32f);
}

void CurveArcLengthMatcher::visualizeMatching(float a_len, float a_off, float b_len, float b_off)
{
	//Get matched subsets of curves
	vector<cv::Point2d> a_subset, b_subset;
	getSubset(source, a_len, a_off, a_subset);
	getSubset(target, b_len, b_off, b_subset);
	ResampleCurve(a_subset, a_subset, RESAMPLE_SIZE, true);
	ResampleCurve(b_subset, b_subset, RESAMPLE_SIZE, true);

	//
	Mat outout(1000, 1000, CV_8UC3, Scalar::all(0));
	drawMatching(outout, a_subset, b_subset);

	imshow("outout", outout);
	//waitKey();
}

void CurveArcLengthMatcher::renderMatching(const string& filename, float a_len, float a_off, float b_len, float b_off)
{
	//Get matched subsets of curves
	vector<cv::Point2d> a_subset, b_subset;
	getSubset(source, a_len, a_off, a_subset);
	getSubset(target, b_len, b_off, b_subset);

	//Normalize to equal length
	ResampleCurve(a_subset, a_subset, RESAMPLE_SIZE, true);
	ResampleCurve(b_subset, b_subset, RESAMPLE_SIZE, true);

	//
	Mat outout(1000, 1000, CV_8UC3, Scalar::all(255));
	drawMatching(outout, a_subset, b_subset);

	//
	//stringstream ss;
	//ss<<this->
	//int fontFace = FONT_HERSHEY_SCRIPT_SIMPLEX;
	//double fontScale = 2;
	//int thickness = 3;
	//int baseline = 0;
	//Size textSize = getTextSize(text, fontFace, fontScale, thickness, &baseline);
	//Point textOrg((outout.cols - textSize.width) / 2, (outout.rows + textSize.height) / 2);
	//putText(outout, text, textOrg, fontFace, fontScale, Scalar::all(255), thickness, 8);

	imwrite(filename.c_str(), outout);
}

void CurveArcLengthMatcher::drawMatching(Mat& outout, vector<cv::Point2d> & a_subset, vector<cv::Point2d> & b_subset)
{
	//Visualize the original and target

	{
		//draw small original
		vector<cv::Point2d> tmp_curve;
		cv::transform(source, tmp_curve, getRotationMatrix2D(Point2f(0, 0), 0, .5));
		Mat tmp_curve_m(tmp_curve); tmp_curve_m += Scalar(5, 0);
		drawOpenCurve(outout, tmp_curve, Scalar(255), 1);

		//draw small matched subset of original
		ConvertCurve(a_subset, tmp_curve);
		cv::transform(tmp_curve, tmp_curve, getRotationMatrix2D(Point2f(0, 0), 0, .5));
		Mat tmp_curve_m1(tmp_curve); tmp_curve_m1 += Scalar(5, 0);
		drawOpenCurve(outout, tmp_curve, Scalar(255, 255), 2);

		//draw small target
		ConvertCurve(target, tmp_curve);
		cv::transform(tmp_curve, tmp_curve, getRotationMatrix2D(Point2f(0, 0), 0, .5));
		Mat tmp_curve_m2(tmp_curve); tmp_curve_m2 += Scalar(outout.cols - 600, 0);
		drawOpenCurve(outout, tmp_curve, Scalar(255, 0, 255), 1);

		//draw small matched subset of target
		ConvertCurve(b_subset, tmp_curve);
		cv::transform(tmp_curve, tmp_curve, getRotationMatrix2D(Point2f(0, 0), 0, .5));
		Mat tmp_curve_m3(tmp_curve); tmp_curve_m3 += Scalar(outout.cols - 600, 0);
		drawOpenCurve(outout, tmp_curve, Scalar(255, 175, 255), 2);

		//draw big target
		drawOpenCurve(outout, target, Scalar(0, 0, 255), 1);
		//draw big matched subset of target
		drawOpenCurve(outout, b_subset, Scalar(0, 255, 255), 2);
	}


	//Prepare the curves for finding the transformation
	Mat trans;
	computeTransform(a_subset, b_subset, trans);
	//cout << trans;

	vector<cv::Point2d> seq_a_trans(source.size());
	cv::transform(source, seq_a_trans, trans);

	vector<cv::Point2d> a_subset_trans(source.size());
	vector<cv::Point2d> tmp_curve;
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



bool CurveArcLengthMatcher::fileExists(const std::string& filename)
{
	struct stat buf;
	if (stat(filename.c_str(), &buf) != -1)
	{
		return true;
	}
	return false;
}
