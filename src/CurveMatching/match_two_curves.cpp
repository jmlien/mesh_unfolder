
#include "CurveMatching/opencvstd.h"
#include "CurveMatching/CurveCSS.h"
#include "CurveMatching/CurveSignature.h"

using namespace cv;

int main(int argc, char** argv)
{
	CurveMatcher matcher;

	///-------------------------------------------------------
	///Read image #1
	///-------------------------------------------------------
	matcher.setSource(argv[1]);
	vector<Point2d> a_p2d;
	ConvertCurve(matcher.getSource(), a_p2d);

	///-------------------------------------------------------
	///Read image #2
	///-------------------------------------------------------
	matcher.setTarget(argv[2]);
	vector<Point2d> b_p2d;
	ConvertCurve(matcher.getTarget(), b_p2d);

	//Compare curves
	int a_len,a_off,b_len,b_off;
	double db_compare_score;
	matcher.CompareCurvesUsingSignatureDB(a_len, a_off, b_len,  b_off, db_compare_score);

	cout << "db_compare_score=" << db_compare_score << endl;

	double wholeRMSE=matcher.ComputeWholeRMSE(a_len, a_off, b_len, b_off);
	cout << "whole model RMSE=" << wholeRMSE << endl;

	matcher.visualizeMatching(a_len, a_off, b_len, b_off);
}
