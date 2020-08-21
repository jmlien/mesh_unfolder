/*
* LaserUnfolding.cpp
*
*  Created on: May 11, 2020
*      Author: jmlien
*/
#include "LaserUnfolding.h"

#include <fstream>
#include <cassert>
#include <ctime>
#include <cfloat>
#include <algorithm>
#include <queue>
#include <sstream>
#include <string>
#include <list>
#include <unordered_set>
using namespace std;

#include <opencv2/core/core.hpp>

#include "UnfoldingProblem.h"
#include "objReader.h"
#include "util/UnfolderHelper.h"
#include "util/DataHelper.h"
#include "util/SVGWriter.h"
#include "util/DisjointSets.h"

using namespace masc::util;
using namespace masc::unfolding::util;


namespace masc {
	namespace unfolding {

		LaserUnfolding::LaserUnfolding(Unfolder* unfolder)
		:BloomingUnfolding(unfolder)
		{
			const Config & mycfg = unfolder->getConfig();
			if(mycfg.unfolding_motion!=Config::Laser_Unfolding){
				cerr << "! Error: LaserUnfolding requires the flag -laser to fold faces in outside-in order" << endl;
				exit(1);
			}

			this->m_no_collision_detection=false;
			this->m_no_laser_visibility=false;
		}

		LaserUnfolding::~LaserUnfolding()
		{

		}

		bool LaserUnfolding::setup(const string& method)
		{
			bool r=BloomingUnfolding::setup(method);

			cout<<"(";
			if (method.find("nocd") != string::npos)
			{
				this->m_no_collision_detection = true; //enable preview only mode
				cout << "\"nocd\" disabled collision detection;";
			}
			else if (method.find("noviz") != string::npos)
			{
				this->m_no_laser_visibility = true; //avoid trimming the model
				cout << "\"noviz\" disabled laser visibility check;";
			}
			else{
				cout << "add \"nocd\" to disable collision detection;";
				cout << "add \"noviz\" to disable laser visibility check;";
			}
			cout<<")"<<endl;


			return true;
		}

		void LaserUnfolding::run()
		{
			cout<<"-------- Optimize for Blooming Unfolding --------"<<endl;
			BloomingUnfolding::run(); //this updates the unfolder
			if(this->m_preview_only || this->m_no_trim) return;

			//optimize for laser forming
			cout<<"---------- Optimize for Laser Unfolding ----------"<<endl;
			float cd_step_size=0.0025; //collision detect step size
			Net laser_net = AStarLaser(cd_step_size);
			m_unfolder->buildFromWeights(laser_net.encode());
			model * mesh = this->m_unfolder->getModel();

			//remove faces that are colliding
			if(!m_no_collision_detection){
				cout<<"- Laser unfolding finding optimal cuts to resolve self-collision"<<endl;

				cout<<"collect_self_intersection 1"<<endl;
				auto overlaps = collect_self_intersection(m_unfolder, laser_net, cd_step_size);
				cout<<"collisions: "; for(int i=0;i<mesh->t_size;i++) if(overlaps[i].empty()==false) cout<<"f["<<i<<"] has "<<overlaps[i].size()<<endl;
				cout<<"collect_self_intersection 2"<<endl;

				BitVector before=laser_net.encode();
				laser_net.optimalcuts(m_unfolder, overlaps);
				BitVector after=laser_net.encode();
				for(int i=0;i<before.size();i++){
					if(before[i]==after[i]) continue;
					cout<<"bit["<<i<<"] changed from "<<before[i]<<" to "<<after[i]<<endl;
					cout<<"edge "<<i<<" has face "<<mesh->edges[i].fid.front()<<" and "<<mesh->edges[i].fid.back()<<endl;
				}

			//	return;

				const vector<uint>& old_ofl=m_unfolder->getOrderedFaceList();
				unordered_map<uint,uint> old2new;

				cout<<"- Laser unfolding building new unfoldings"<<endl;
				m_unfolder = build_unfolder_from_net(m_unfolder, laser_net);
				model * subm=m_unfolder->getModel();
				for(int i=0;i<subm->t_size;i++)
					old2new[ subm->tris[i].source_fid ]=i;
				vector<uint> new_ofl;
				for(uint fid : old_ofl){
					if(old2new.find(fid)!=old2new.end()) new_ofl.push_back(old2new[fid]);
				}
				assert(new_ofl.size()==subm->t_size);
				m_unfolder->setOrderedFaceList(new_ofl);

				stringstream ss;
				ss << "./" << mesh->name << "_laser_" << subm->t_size << ".obj";
				subm->saveObj(ss.str());
				cout<<"- Laser unfolding saved the trimmed model to "<<ss.str()<<endl;


				cout<<"collect_self_intersection 3"<<endl;
				overlaps = collect_self_intersection(m_unfolder, laser_net, cd_step_size);
				cout<<"collisions: "; for(int i=0;i<overlaps.size();i++){
					if(overlaps[i].empty()==false){
						cout<<"f["<<i<<", aka "<<subm->tris[i].source_fid<<"] has "<<overlaps[i].size()<<endl;
					}
				}
				cout<<"collect_self_intersection 4"<<endl;

			}

			//check for laser visibility


		}// END LaserUnfolding::run()


		LaserUnfolding::Net LaserUnfolding::AStarLaser(float cd_step)
		{
			Net start=getNet(m_unfolder);
			AStar_Helper_Laser laser(start, m_unfolder, cd_step, 1, 1, 100); //stop when there is 0 overlap and try it until level_3
			return AStar(laser);
		}


		//collect all intersections and place them inside the unfolder
		vector<set<uint>> LaserUnfolding::collect_self_intersection
		(Unfolder * unfolder, const LaserUnfolding::Net& net, float cd_step)
		{
				//return the number of colliding face during the folding process

				//cout<<"- Evaluate folding motion....."<<endl;
				unordered_set<int> colliding;
				model * mesh=unfolder->getModel();
				vector<set<uint>> collected_overlaps;
				collected_overlaps.resize(mesh->t_size);

				for(float i=1.0;i>=0;i-=cd_step)
				{
					unfolder->unfoldTo(i);
					int count=unfolder->checkCollision();
					if(count>0){
						const vector<set<uint>>& overlaps = unfolder->getOverlppingFacePairs();
						for(int j=0;j<mesh->t_size;j++){
							if(overlaps[j].empty()) continue;
							for(int o : overlaps[j]) collected_overlaps[j].insert(o);
						}//end for j
					}//end if
				}//end for i

				return collected_overlaps;
		}


		//
		//
		// A* sub-modules
		//
		//

		//-------------------
		//in this helper, we only flip the edges so there without affecting the keep faces

		float LaserUnfolding::AStar_Helper_Laser::dist
		(LaserUnfolding::Net& net, const vector<float>& ecosts)
		{
				model * m=unfolder->getModel();
				if (net.G() == FLT_MAX || g_recompute_dist) {
					float flips = 0;
					for (int i = 0; i < net.encode().size(); i++){
						if (this->root.encode()[i] != net.encode()[i]) {
							//flips++;
							const edge& e=m->edges[i];
							flips+=e.length;
							//cout<<"Ahhhh..."<<endl;
						}
					}
					net.G()=flips;
				}

				return net.G();
		}

		float LaserUnfolding::AStar_Helper_Laser::heuristics
		(LaserUnfolding::Net& net)
		{

				if (net.H() == FLT_MAX) {
					float h=0;
					//number of overlaps
					int n=unfolder->buildFromWeights(net.encode());

					if(n>0) h=FLT_MAX/2; //no overlapping allowed
					else h=evaluate_folding_motion(unfolder);
					net.H()=h;
				}

				return net.H();
		}


		//return the number of colliding face during the folding process
		int LaserUnfolding::AStar_Helper_Laser::evaluate_folding_motion(Unfolder * unfolder)
		{
					//cout<<"- Evaluate folding motion....."<<endl;
					unordered_set<int> colliding;
					model * mesh=unfolder->getModel();

					for(float i=1.0;i>=0;i-=m_cd_step)
					{
						unfolder->unfoldTo(i);
						if (i >= 0.990001){
							int count=unfolder->checkOverlaps();
							//cout<<count<<" collision found... at step="<<i<<endl;
							if(count>0){
								const vector<set<uint>>& overlaps = unfolder->getOverlppingFacePairs();
								for(int i=0;i<mesh->t_size;i++){
									if(overlaps[i].empty()) continue;
									colliding.insert(i);
									//cout<<"face["<<i<<"] overlaps with";
									//for(int o : overlaps[i]) cout<<o<<";";
									//cout<<endl;
								}
							}
						}
						else{
								int count=unfolder->checkCollision();
								//cout<<count<<" collision found... at step="<<i<<endl;
								if(count>0){
									for(int j=0;j<mesh->t_size;j++){
										triangle& t=mesh->tris[j];
										if(!t.overlapped) continue;
										colliding.insert(j);
										//cout<<"face["<<j<<"] overlapped"<<endl;
									}//end for j
								}
						}
					}//end for i

				/*	cout<<"\t- Colliding faces: ";
					for(int i : colliding) cout<<i<<";";
					cout<<endl;*/

					return colliding.size();
		}



	}
}
