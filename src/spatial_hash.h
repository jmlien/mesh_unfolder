#pragma once


#include <unordered_map>
#include <list>
#include <vector>
#include <cassert>
using namespace std;

#include <limits.h>
#include <Point.h>
#include <Vector.h>
using namespace mathtool;

typedef unsigned int uint;

struct spatial_hash
{
	spatial_hash(const mathtool::Point3d& origin, double w = 1, double h = 1, double d = 1)
	{
		hash_O = origin;
		hash_w = w;
		hash_h = h;
		hash_d = d; 
	}

	uint add_point(const mathtool::Point3d& pos);  //add a point to hast
	uint find_point(const mathtool::Point3d& pos); //query/find the id of a given point
	const vector<mathtool::Point3d> & get_all_points(){ return all_points; }

private:

	//hash_func
	uint spatial_hash_func(const mathtool::Point3d& pos);
	unordered_map<uint, list<uint> > hash;

	//these are used to hash points
	mathtool::Point3d hash_O;
	double  hash_w, hash_h, hash_d; //width (X), height (Y), depth (Z)

	vector<mathtool::Point3d> all_points;
};

