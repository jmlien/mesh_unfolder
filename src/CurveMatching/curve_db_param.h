/*
 * curve_db_param.h
 *
 *  Created on: Oct 1, 2015
 *      Author: zxi
 */

#ifndef CURVE_DB_PARAM_H_
#define CURVE_DB_PARAM_H_


//-------------------------------------------------------------------------------------------------------------------------------
//
// Varibles and Functions for controlling the structure/size of the curve database
// Note these variables are all global state varibles.
//
template<typename MATCHER>
struct CURVE_DB_PARAM
{
	static int CurveResampleSize;
	static typename MATCHER::position_type SmallestCurveSegmentSize;
	static typename MATCHER::position_type LongestCurveSegmentSize;
	static typename MATCHER::position_type OffsetStepSize;
	static int CurvatureFilterSize;

	static void initializeMatcher(MATCHER & matcher)
	{
		matcher.setCurveResampleSize(CURVE_DB_PARAM<MATCHER>::CurveResampleSize);
		matcher.setSmallestCurveSegmentSize(CURVE_DB_PARAM<MATCHER>::SmallestCurveSegmentSize);
		matcher.setLongestCurveSegmentSize(CURVE_DB_PARAM<MATCHER>::LongestCurveSegmentSize);
		matcher.setOffsetStepSize(CURVE_DB_PARAM<MATCHER>::OffsetStepSize);
		matcher.setCurvatureFilterSize(CURVE_DB_PARAM<MATCHER>::CurvatureFilterSize);
	}
};

template<typename MATCHER> int CURVE_DB_PARAM<MATCHER>::CurveResampleSize = 100;
template<typename MATCHER> typename MATCHER::position_type CURVE_DB_PARAM<MATCHER>::SmallestCurveSegmentSize = 50;
template<typename MATCHER> typename MATCHER::position_type CURVE_DB_PARAM<MATCHER>::LongestCurveSegmentSize = 75;
template<typename MATCHER> typename MATCHER::position_type CURVE_DB_PARAM<MATCHER>::OffsetStepSize = 2;
template<typename MATCHER> int CURVE_DB_PARAM<MATCHER>::CurvatureFilterSize = 1000;



#endif /* CURVE_DB_PARAM_H_ */
