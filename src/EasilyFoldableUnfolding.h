/*
* EasilyFoldableUnfolding.h
*
* Create Linearly or Uniformly Foldable Net using A* strategies
*
*  Created on: Aug 21, 2020
*      Author: jmlien
*/

#pragma once

#include <string>
#include <vector>

#include "BloomingUnfolding.h"
#define LAZY

namespace masc {
	namespace unfolding {

		using namespace std;

		class EasilyFoldableUnfolding : public BloomingUnfolding
		{
		public:

			EasilyFoldableUnfolding(Unfolder* unfolder);
			virtual ~EasilyFoldableUnfolding();

			bool setup(const string& method) override;

			//////////////////////////////////////////
			// run n unfoldings and cluster
			//////////////////////////////////////////

			void run() override;

			//////////////////////////////////////////
			// access
			//////////////////////////////////////////

		protected:

			Net AStarFoldable(float cd_step);

			struct AStar_Helper_Foldable : public AStar_Helper_Mixed
			{

				AStar_Helper_Foldable(const Net& net, Unfolder * uf, float step, float t=0, int ml=3, int mf=1000)
				:AStar_Helper_Mixed(net,uf,t,ml,mf), m_cd_step(step){}

				list<Net> neighbors(Net& net) override;
				float dist(Net& net, const vector<float>& ecosts) override;
				float heuristics(Net& net) override;
				int evaluate_folding_motion(Unfolder * unfolder);
				void find_path(Net& net, list<uint>& path);
				void find_path(Net& net, uint f1, uint f2, list<uint>& path);

				//collect all intersecting pairs during the entire folding process
				vector<set<uint>> collect_self_intersection(const Net& net);

				float m_cd_step;
			};

#ifdef LAZY
			//optimize the net with A*
			Net AStar(AStar_Helper & help) override;
#endif

			bool m_no_collision_detection;
			bool m_no_laser_visibility;

		};


	}
}
