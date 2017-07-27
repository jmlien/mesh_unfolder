#include "CurveMatching/img2ply.h"

#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>

//using namespace cv;
using namespace std;

inline void contour2ply(vector<cv::Point>& contour, c_ply& ply)
{
	int psize = contour.size();
	ply.beginPoly();
	for (int j = 0; j < psize; j++)
		ply.addVertex(contour[j].x, contour[j].y);
	ply.endPoly();
}

bool img2ply(cv::Mat src, c_polygon& polygon)
{
	  cv::Mat src_gray;
	  int thresh = 50;

		//negate the color
		bitwise_not(src, src_gray);

		/// Convert image to gray and blur it
		//cv::cvtColor(src, src_gray, CV_BGR2GRAY);
		//blur( src_gray, src_gray, Size(3,3) );



		/// Detect edges using Threshold
		cv::Mat threshold_output;
		cv::threshold(src_gray, threshold_output, thresh, 255, cv::THRESH_BINARY);


		//cv::namedWindow("img", CV_WINDOW_AUTOSIZE);
		//cv::imshow("img", threshold_output);
		//cv::waitKey(0);                                          // Wait for a keystroke in the window


		/// Find contours

		vector<vector<cv::Point> > contours;
		vector<cv::Vec4i> hierarchy;

		cv::findContours(threshold_output, contours, hierarchy, CV_RETR_TREE, CV_CHAIN_APPROX_SIMPLE);
		//cv::findContours(threshold_output, contours, CV_RETR_EXTERNAL, CV_CHAIN_APPROX_SIMPLE);

		int csize = contours.size();

		//find the largest external boundary
		int outid = -1;
		int max_size = -1;

		//For each i - th contour contours[i], the elements hierarchy[i][0], hiearchy[i][1], hiearchy[i][2], and hiearchy[i][3] are set to 0 - based indices in
		//contours of the next and previous contours at the same hierarchical level, the first child contour and the parent contour, respectively.If for the
		//contour i there are no next, previous, parent, or nested contours, the corresponding elements of hierarchy[i] will be negative.

		for (int i = 0; i < csize; i++)
		{
			//cout << "contour[" << i << "] size = " << contours[i].size() << " has hie[0]=" << hierarchy[i][0] << " hie[1]=" << hierarchy[i][1] << " hie[2]=" << hierarchy[i][2] << " hie[3]=" << hierarchy[i][3] << endl;

			if (hierarchy[i][3] == -1) //external boundary
			{
				int csize = contours[i].size();
				if (csize>max_size)
				{
					max_size = csize;
					outid = i;
				}
			}
		}//end for i

		if (outid < 0) return false;

		//------------------------------------------------------------------
		//create the external boundary
		c_ply ply(c_ply::POUT);
		contour2ply(contours[outid],ply);
		polygon.push_back(ply);

		//create the children of the external boundary
		int childid = hierarchy[outid][2];
		while (childid != -1)
		{
			c_ply ply(c_ply::PIN);
			contour2ply(contours[childid], ply);
			polygon.push_back(ply);
			childid = hierarchy[childid][0]; //move to the previous child
		}

		return true;
}

//given a black/whit image with size (widthxheight) return a ploygon representing the contour
bool img2ply(unsigned char * img, unsigned int width, unsigned int height, c_polygon& polygon)
{
	/// Load source image and convert it to gray
	cv::Mat src = cv::Mat(height, width, CV_8UC1, img); //imread(argv[1], 1);
	return img2ply(src,polygon);
}


//Note: this function extracts all contours except those with size smaller than threshold
bool img2ply_complete(unsigned char * img, unsigned int width, unsigned int height, vector<c_polygon>& polygons, float threshold)
{
	cv::Mat src_gray;
	int thresh = 50;

	/// Load source image and convert it to gray
	cv::Mat src = cv::Mat(height, width, CV_8UC1, img); //imread(argv[1], 1);

	//negate the color
	bitwise_not(src, src_gray);

	/// Detect edges using Threshold
	cv::Mat threshold_output;
	cv::threshold(src_gray, threshold_output, thresh, 255, cv::THRESH_BINARY);

	return img2ply_complete(threshold_output, polygons, threshold);
}


//Note: this function extracts all contours except those with size smaller than threshold
bool img2ply_complete(cv::Mat& img, vector<c_polygon>& polygons, float threshold)
{

	/// Find contours

	vector<vector<cv::Point> > contours;
	vector<cv::Vec4i> hierarchy;

	cv::findContours(img, contours, hierarchy, CV_RETR_TREE, CV_CHAIN_APPROX_SIMPLE);

	int csize = contours.size();

	//For each i - th contour contours[i], the elements hierarchy[i][0], hiearchy[i][1], hiearchy[i][2], and hiearchy[i][3] are set to 0 - based indices in
	//contours of the next and previous contours at the same hierarchical level, the first child contour and the parent contour, respectively.If for the
	//contour i there are no next, previous, parent, or nested contours, the corresponding elements of hierarchy[i] will be negative.

	for (int i = 0; i < csize; i++)
	{
		if (hierarchy[i][3] == -1) //external boundary
		{
			int csize = contours[i].size();
			//if (csize <= threshold) continue; //this external boundary is too small... ignore

			c_polygon polygon;

			//------------------------------------------------------------------
			//create the external boundary
			c_ply ply(c_ply::POUT);
			contour2ply(contours[i], ply);
			polygon.push_back(ply);

			//create the hole boundaries of the external boundary
			int childid = hierarchy[i][2];

			while (childid != -1)
			{
				c_ply ply(c_ply::PIN);
				if (contours[childid].size() >= threshold) //make sure that the hole boundary is big enough6
				{
					contour2ply(contours[childid], ply);
					polygon.push_back(ply);
				}
				childid = hierarchy[childid][0]; //move to the previous child
			}
			//------------------------------------------------------------------

			if (polygon.getArea()>1e-10 && polygon.getArcLength()>threshold ) //at least this polygon needs to have non-zero area
				polygons.push_back(polygon); //done, store the result
		}
	}//end for i

	return true;
}
