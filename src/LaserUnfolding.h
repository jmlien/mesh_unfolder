/*
* LaserUnfolding.h
*
*  Created on: May 11, 2020
*      Author: jmlien
*/

#pragma once

#include <string>
#include <vector>

#include "BloomingUnfolding.h"

namespace masc {
	namespace unfolding {

		using namespace std;

		class LaserUnfolding : public BloomingUnfolding
		{
		public:

			LaserUnfolding(Unfolder* unfolder);
			virtual ~LaserUnfolding();

			bool setup(const string& method) override;

			//////////////////////////////////////////
			// run n unfoldings and cluster
			//////////////////////////////////////////

			void run() override;

			//////////////////////////////////////////
			// access
			//////////////////////////////////////////

		protected:

			//collect all intersecting pairs during the entire folding process
			vector<set<uint>> collect_self_intersection
			(Unfolder * unfolder, const Net& net, float cd_step);

			Net AStarLaser(float cd_step);

			struct AStar_Helper_Laser : public AStar_Helper_Mixed
			{

				AStar_Helper_Laser(const Net& net, Unfolder * uf, float step, float t=0, int ml=3, int mf=1000)
				:AStar_Helper_Mixed(net,uf,t,ml,mf), m_cd_step(step){}
				float dist(Net& net, const vector<float>& ecosts) override;
				float heuristics(Net& net) override;
				int evaluate_folding_motion(Unfolder * unfolder);
				float m_cd_step;
			};

			bool m_no_collision_detection;
			bool m_no_laser_visibility;

		};


	}
}
