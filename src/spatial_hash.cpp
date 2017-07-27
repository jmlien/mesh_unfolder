#include "spatial_hash.h"
#include "intersection.h"

uint spatial_hash::spatial_hash_func(const Point3d& pos)
{
	uint x = ((pos[0] - hash_O[0])*1.0 / hash_w) * 100;
	uint y = ((pos[1] - hash_O[1])*1.0 / hash_h) * 100;
	uint z = ((pos[2] - hash_O[2])*1.0 / hash_d) * 100;

	uint offset=sizeof(uint) / 2;
	return (((x << offset) | y) << offset) | z;
}

bool AlmostEqual3(const Point3d& a, const Point3d& b, double smallnumber)
{
	return (fabs(a[0] - b[0])<smallnumber &&fabs(a[1] - b[1])<smallnumber && fabs(a[2] - b[2])<smallnumber);
}

uint spatial_hash::add_point(const Point3d& pos)  //add a point to hast
{
	uint key = spatial_hash_func(pos);
	hash_map::iterator it = hash.find(key);

	if (it != hash.end())
	{
		list<uint>& collision = it->second;

		for (list<uint>::iterator i = collision.begin(); i != collision.end(); i++)
		{
			uint id = *i;
			bool almosteq = AlmostEqual3(all_points[id], pos, 0.01);
			if (almosteq){ return id; }
		}
	}
	
	//not in hash...
	//add to node vector
	uint id = all_points.size();
	all_points.push_back(pos);

	//add to hash
	if (it != hash.end()) it->second.push_back(id);
	else{
		list<uint> tmp;
		tmp.push_back(id);
		hash[key] = tmp;
	}

	return id;
}

uint spatial_hash::find_point(const Point3d& pos) //query/find the id of a given point
{
	uint key = spatial_hash_func(pos);
	hash_map::iterator it = hash.find(key);

	if (it != hash.end())
	{
		list<uint>& collision = it->second;

		for (list<uint>::iterator i = collision.begin(); i != collision.end(); i++)
		{
			uint id = *i;
			bool almosteq = AlmostEqual(all_points[id].get(), pos.get());
			if (almosteq){ return id; }
		}
	}
	
	return UINT_MAX;
}