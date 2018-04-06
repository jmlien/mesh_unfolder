/*
 *  Software License Agreement (BSD License)
 *
 *  Copyright (c) 2012, Willow Garage, Inc.
 *  All rights reserved.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions
 *  are met:
 *
 *   * Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *   * Redistributions in binary form must reproduce the above
 *     copyright notice, this list of conditions and the following
 *     disclaimer in the documentation and/or other materials provided
 *     with the distribution.
 *   * Neither the name of Willow Garage, Inc. nor the names of its
 *     contributors may be used to endorse or promote products derived
 *     from this software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 *  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 *  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *  COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 *  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 *  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 *  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 *  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 *  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 *  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 *
 *  File:    spectralClustering.hpp
 *  Author:  Hilton Bristow
 *  Created: Aug 24, 2012
 */

#ifndef SPECTRALCLUSTERING_HPP_
#define SPECTRALCLUSTERING_HPP_

#include <limits>
#include <cstdlib>
#include <algorithm>
#include <iostream>
#include <opencv2/core/core.hpp>
#include <opencv2/core/core_c.h>
using namespace std;

#ifdef __cplusplus

/*! \namespace cv
Namespace where all the C++ OpenCV functionality resides
*/
namespace cv
{

	// ---------------------------------------------------------------------------
	// DISTANCE FUNCTIONS
	// ---------------------------------------------------------------------------

	/*! /class DistanceFunction
	 *  /brief A functor interface which defines a distance or similarity function
	 *
	 *	This class provides a standard functor interface for distance and
	 *	similarity functions (which are both considered to be members of
	 *	comparison classes). Subclasses of this function can then be used in
	 *	any function that requires as input a DistanceFunction.
	 *
	 *	The DistanceFunction has a single public data member, D, the dimensionality
	 *	of the data, and an operator() which performs the distance computation of
	 *	two input data arrays.
	 *
	 *	Subclasses are also encouraged to provide a static method with the signature
	 *	T compute(T const * const p1, T const * const p1, const uint8_t D, ...), which
	 *	constructs an instance and calls the operator() for once-off distance
	 *	computations.
	 *
	 */
	//template<typename T>
	//class DistanceFunction {
	//public:
	//	/// the dimensionality
	//	const uint8_t D;
	//	/// constructor for initialising the data dimensionality
	//	DistanceFunction(const uint8_t _D) : D(_D) {}
	//	/// virtual destructor
	//	virtual ~DistanceFunction() {};

	//	/*! /brief distance operator
	//	 *
	//	 * When called, computes the distance between points p1 and p2
	//	 * using the implemented distance measure. This assumes that
	//	 * p1[D-1] and p2[D-1] point to valid memory addresses for speed
	//	 * since this function will often be called in a loop.
	//	 *
	//	 * @param p1 the first point
	//	 * @param p2 the second point
	//	 * @return the implemented measure of distance or similarity between
	//	 * p1 and p2
	//	 */
	//	virtual T operator() (T const * const p1, T const * const p2) const = 0;
	//};

	///*! /class EuclideanEpsilonSimilarity
	// *  /brief Euclidean similarity function for spectral clustering
	// *
	// *  This class implements the DistanceFunction interface using Euclidean
	// *  distance. Since spectral clustering requires similarity functions,
	// *  this uses the binary epsilon method, such that:
	// *
	// *  sqrt(sum( (p2 - p1)^2 )) < eps ? 1 : 0;
	// */
	//template<typename T>
	//class EuclideanEpsilonSimilarity : public DistanceFunction<T> {
	//public:
	//	const T eps;
	//	EuclideanEpsilonSimilarity(const uint8_t _D, const T _eps) : DistanceFunction<T>(_D), eps(_eps) {}
	//	T operator() (T const * const p1, T const * const p2) const {
	//		T dist = 0;
	//		for (uint8_t n = 0; n < this->D; ++n) dist += pow(p2[n] - p1[n], 2);
	//		return sqrt(dist) < eps ? 1 : 0;
	//	}
	//	static T compute(T const * const p1, T const * const p2, const uint8_t D, const T eps) {
	//		return EuclideanEpsilonSimilarity(D, eps).operator()(p1, p2);
	//	}
	//	static EuclideanEpsilonSimilarity TwoD(T _eps) { return EuclideanEpsilonSimilarity(2, _eps); }
	//	static EuclideanEpsilonSimilarity ThreeD(T _eps) { return EuclideanEpsilonSimilarity(3, _eps); }

	//};




	//// ---------------------------------------------------------------------------
	//// SPECTRAL CLUSTERING
	//// ---------------------------------------------------------------------------

	///*! /brief compute the pariwise similarity matrix for a set of points
	// *
	// * Given a set of datapoints, data, where observations are rows and variables
	// * are columns, and a similarity measure implementing the DistanceFunction
	// * interface, this returns the adjacency matrix for the data.
	// *
	// * This method is intended to construct input for the spectralClustering
	// * method, and thus asserts that the distance is symmetric.
	// *
	// * @param data the input datapoints
	// * @param adjacency the output adjacency matrix
	// * @param distance the distance function
	// */
	//template<typename T>
	//void adjacencyMatrix(const cv::Mat_<T>& data, cv::Mat_<T>& adjacency, DistanceFunction<T>& distance) {

	//	// assert that the distance function is symmetric by random sampling
	//	unsigned int D = data.size().width;
	//	unsigned int N = data.size().height;
	//	CV_Assert(distance.D == D);
	//	for (unsigned int n = 0; n < std::min((unsigned int)5, N); ++n) {
	//		unsigned int i = rand() % N;
	//		unsigned int j = rand() % N;
	//		CV_Assert(distance(data[i], data[j]) == distance(data[j], data[i]));
	//	}

	//	// construct the upper triangular part of the adjacency matrix
	//	adjacency.create(N, N);
	//	for (unsigned int i = 0; i < N; ++i) {
	//		for (unsigned int j = i; j < N; ++j) {
	//			adjacency(i, j) = distance(data[i], data[j]);
	//		}
	//	}
	//	// replicate to produce symmetry
	//	completeSymm(adjacency);
	//}

	// the normalization schemes available
	enum {
		LAPLACIAN_UNNORMALIZED = 0,
		LAPLACIAN_SHI_MALIK = 1,
		LAPLACIAN_NG_WEISS = 2
	};

	/*! /brief compute the graph laplacian from an adjacency matrix
	 *
	 * Given an adjacency matrix and a normalization method, compute the
	 * graph laplacian. The normalization method can be any one of:
	 * 	- LAPLACIAN_UNNORMALIZED
	 * 		- no normalization applied
	 * 	- LAPLACIAN_SHI_MALIK
	 * 		- normalization according to Shi and Malik, "Normalized Cuts
	 * 		  and Image Segmentation", PAMI 2000
	 * 	- LAPLACIAN_NG_WEISS
	 * 		- normalization according to Ng, Jordan and Weiss,
	 * 		  "On Spectral Clustering: Analysis and an Algorithm", NIPS 2002
	 *
	 * @param W the adjacency matrix
	 * @param L the graph laplacian
	 * @param normalization the normalization method
	 */
	template<typename T>
	void graphLaplacian(const cv::Mat_<T>& W, cv::Mat_<T>& L, const int normalization) {

		// compute the degree matrix
		unsigned int N = W.rows;
		cv::Mat_<T> D;
		cv::reduce(W, D, 0, CV_REDUCE_SUM);

		// compute the laplacian
		switch (normalization) {
		case LAPLACIAN_UNNORMALIZED: {
			L = cv::Mat::diag(D) - W;
			break;
		} case LAPLACIAN_NG_WEISS: {
			cv::Mat_<T> I = cv::Mat_<T>::eye(N, N);
			cv::pow(D, -0.5, D);
			D = cv::Mat::diag(D);
			L = I - (D * W * D);
			break;
		} case LAPLACIAN_SHI_MALIK: {
			cv::Mat_<T> I = cv::Mat_<T>::eye(N, N);
			L = I - cv::Mat::diag(1 / D) * W;
			break;
		} default:
			CV_Error(CV_StsBadArg, "Unknown graph laplacian normalization method supplied");
			break;
		}
	}

	template<typename T> string serialization_filename(const cv::Mat_<T>& A)
	{
		stringstream ss;
		ss.precision(3);
		auto w = A.rows;
		auto h = A.cols;

		double min, max;
		cv::Point min_loc, max_loc;
		cv::minMaxLoc(A, &min, &max, &min_loc, &max_loc);

		ss << ".spectral_" << w << "_" << h << "_"
			               << min << "_" << max << "_"
						   << min_loc.x << "_" << min_loc.y << "_"
						   << max_loc.x << max_loc.y << ".txt";
		return ss.str();
	}

	inline bool file_exist(const string& filename) {
		if (FILE *file = fopen(filename.c_str(), "r")) {
			fclose(file);
			return true;
		}
		else
			return false;
	}

	template<typename T> bool save_eigin_matrices(const cv::Mat_<T>& A, const cv::Mat_<T>& L, const cv::Mat_<T>& V, const cv::Mat_<T>& D)
	{
		string fname = serialization_filename(A);
		cerr << "- Save to " << fname << endl;

		cv::FileStorage fs(fname, cv::FileStorage::WRITE);
		fs << "A" << A;
		fs << "L" << L;
		fs << "V" << V;
		fs << "D" << D;
		fs.release();
		return true;
	}

	template<typename T> bool load_eigin_matrices(const cv::Mat_<T>& A, cv::Mat_<T>& L, cv::Mat_<T>& V, cv::Mat_<T>& D)
	{
		string filename = serialization_filename(A);
		if (file_exist(filename) == false)
		{
			cerr << "! Warning: Failed loading " << filename << endl;
			return false;
		}
		cv::FileStorage fs(filename, cv::FileStorage::READ);
		cv::Mat_<T> A2;

		fs["A"] >> A2;
		fs["L"] >> L;
		fs["V"] >> V;
		fs["D"] >> D;
		fs.release();

		cv::Mat diff = (A != A2); //make sure A and A2 ar identical
		return cv::countNonZero(diff) == 0;
	}

	/*! @brief perform spectral clustering
	 *
	 * Given the supplied adjacency matrix, the number of clusters and a graph laplacian
	 * normalization method, compute the clusters on the smallest k eigenvectors of
	 * the graph laplacian.
	 *
	 * The adjacency matrix can be computed from raw data using the adjacencyMatrix()
	 * method. The current similarity measure implemented is binary epsilon similarity
	 * using euclidean distance (EuclideanEpsilonDistance). If you wish to use a different
	 * similarity measure, simply provide your own class derived from DistanceFunction.
	 *
	 * @param A the adjacency matrix (N x N)
	 * @param idx the output matrix of cluster membership (N x 1)
	 * @param K the number of desired clusters for the kmeans step
	 * @param normalize the graph laplacian normalization method
	 */
	template<typename T>
	void spectralClustering(const cv::Mat_<T>& A, cv::Mat& idx, unsigned int K, int normalize = LAPLACIAN_UNNORMALIZED)
	{

		cv::Mat_<T> L, V, D;
		if (load_eigin_matrices(A,L,V,D)==false)
		{
			cout << "compute the Laplacian of the adjacency matrix" << endl;
			graphLaplacian(A, L, normalize);

			// calculate the eigenvalues and eigenvectors
			cout << "calculate the eigenvalues and eigenvectors" << endl;
			cv::eigen(L, D, V);

			// keep only the k-principal components
			cout << "keep only the k-principal components" << endl;
			unsigned int N = V.rows;
			V = V.rowRange(N - K - 1, N - 1);

			save_eigin_matrices(A, L, V, D);
		}

		// if using the Ng, Jordan and Weiss laplacian, normalize the eigenvectors
		if (normalize == LAPLACIAN_NG_WEISS) {
			cout << "normalize the eigenvectors" << endl;
			cv::Mat_<T> norm;
			cv::pow(V, 2, norm);
			cv::reduce(norm, norm, 0, CV_REDUCE_SUM);
			for (unsigned int k = 0; k < K; ++k) V.row(k) = V.row(k) / norm;
		}

		// cluster the results
		cv::Mat_<T> centers;
		cv::Mat_<float> Vf;
		cv::TermCriteria term(1, 100, 1e-9);
		V = V.t();
		V.convertTo(Vf, CV_32F);
		kmeans(Vf, K, idx, term, 5, cv::KMEANS_PP_CENTERS, centers);
	}
}//end namespace cv

#endif /* __cplusplus */

#endif /* SPECTRALCLUSTERING_HPP_ */
