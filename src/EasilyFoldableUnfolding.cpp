/*
* EasilyFoldableUnfolding.cpp
*
* Create Linearly or Uniformly Foldable Net
*
*  Created on: Aug 21, 2020
*      Author: jmlien
*/
#include "EasilyFoldableUnfolding.h"

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

		list<int> neighbor_edges2(const BitVector& code, model * m, int eid);

		EasilyFoldableUnfolding::EasilyFoldableUnfolding(Unfolder* unfolder)
		:BloomingUnfolding(unfolder)
		{
			const Config & mycfg = unfolder->getConfig();

			this->m_no_collision_detection=false;
			this->m_no_laser_visibility=false;
		}

		EasilyFoldableUnfolding::~EasilyFoldableUnfolding()
		{

		}

		bool EasilyFoldableUnfolding::setup(const string& method)
		{
			bool r=BloomingUnfolding::setup(method);

			// cout<<"(";
			// if (method.find("nocd") != string::npos)
			// {
			// 	this->m_no_collision_detection = true; //enable preview only mode
			// 	cout << "\"nocd\" disabled collision detection;";
			// }
			// else if (method.find("noviz") != string::npos)
			// {
			// 	this->m_no_laser_visibility = true; //avoid trimming the model
			// 	cout << "\"noviz\" disabled laser visibility check;";
			// }
			// else{
			// 	cout << "add \"nocd\" to disable collision detection;";
			// 	cout << "add \"noviz\" to disable laser visibility check;";
			// }
			// cout<<")"<<endl;


			return true;
		}

		void EasilyFoldableUnfolding::run()
		{
			//cout<<"-------- Optimize for Blooming Unfolding --------"<<endl;
			//BloomingUnfolding::run(); //this updates the unfolder
			//if(this->m_preview_only || this->m_no_trim) return;
			const auto config = this->m_unfolder->getConfig();
			model * mesh = this->m_unfolder->getModel();

			m_start_time = clock();

			//unfolding using blooming unfokding
			blooming_unfold(m_unfolder);

			//compute the face types from a given vector
			//in this caes the normal of the base face
			//compute_face_types(m_unfolder);

			if (this->m_preview_only) return; //preview the unfolding without optimization

			//optimize for Easily Foldable Unfolding
			cout<<"---------- Optimize for Easily Foldable Unfolding ----------"<<endl;
			float cd_step_size=0.4; //collision detect step size
			Net folable_net = AStarFoldable(cd_step_size);
			m_unfolder->buildFromWeights(folable_net.encode());

		}// END EasilyFoldableUnfolding::run()


		EasilyFoldableUnfolding::Net EasilyFoldableUnfolding::AStarFoldable(float cd_step)
		{
			Net net=getNet(m_unfolder);
#ifdef LAZY
			AStar_Helper_Foldable helper(net, m_unfolder, cd_step, 0, 6);
#else
			AStar_Helper_Foldable helper(net, m_unfolder, cd_step, 0, 1);
#endif

			return this->AStar(helper);
		}


#ifdef LAZY
		//optimize the net with Lazy A*
		EasilyFoldableUnfolding::Net EasilyFoldableUnfolding::AStar(AStar_Helper & help)
		{
			//build the first net from unfolder
			Unfolder * unfolder=help.unfolder;
			model * m = unfolder->getModel();
			const int fbase = unfolder->getBaseFaceID(); //base face of a net, same of all nets

			//compute additional data to help guiding A* search
			unordered_set<BitVector,hash_BitVector> visted;
			vector<float> ecosts; //cost of flipping a root edge
			help.root.compute_edge_cost(m, ecosts);

			help.heuristics(help.root);
			help.dist(help.root,ecosts);
			visted.insert(help.root.encode());


			Net best_net = help.root;
			{//get best base faces with min H value
				Config mycfg = unfolder->getConfig();
				int best_fid=fbase;
				mycfg.random_baseface=false;
				//int base=unfolder->getBaseFaceID();
				for(int fid=0;fid<m->t_size;fid++){
					if(fbase==fid) continue;
					mycfg.baseface=fid;
					Unfolder tmp(m,mycfg);
					blooming_unfold(&tmp);
					Net n=getNet(&tmp);
					help.heuristics(n);
					help.dist(n,ecosts);
					if(best_net.H()>n.H()){ best_net=n; best_fid=fid; }
					visted.insert(n.encode());
				}
				cout<<"- Best base face = "<<best_fid<<endl;
			}

		  help.root = best_net;
			vector<Net> best_nets(help.max_level+1,best_net);

			//start A* seach from the start net
			vector<vector<Net>> opens(help.max_level, vector<Net>());
			opens[0].push_back(help.root);
			make_heap(opens[0].begin(), opens[0].end());

			vector<int> failed_attempts(help.max_level+1,0);

			cout << "- Start Lazy A* optimization with "<< help.max_level<< " levels"<<endl;
			int best_level_max=1;
			int current_level=0;

			while (true) {

				//cout<<"current_level="<<current_level<<endl;
				vector<Net>& open=opens[current_level];
				if(open.empty()){
					current_level=(current_level+1)%help.max_level;
					continue;
				}

				Net net = open.front();
				pop_heap(open.begin(), open.end());
				open.pop_back();

				hash_BitVector bvhash;
				//cout<<"Current level="<<net.Level()<<" G="<<net.G()<<" H="<<net.H()<<" open size="<< open.size()<<endl; // << " hash="<<bvhash(net.encode())<<endl;

				//determine if we are making progress
				if ( (help.heuristics(net) < help.heuristics(best_nets[net.Level()])) ||
				( (help.heuristics(net)== help.heuristics(best_nets[net.Level()])) && help.dist(net,ecosts)<help.dist(best_nets[net.Level()],ecosts)) )
				{
						cout << "- Best net overlaps=" << help.heuristics(net) <<", level="<<net.Level()
							<< " dist=" << help.dist(net,ecosts) << ", open size=" << open.size() << endl;
						//reset only when the overall score is descreased
						//if (help.heuristics(net)+help.dist(net,ecosts) <= help.heuristics(best_net)+help.dist(best_net,ecosts))
						failed_attempts[net.Level()] = 0; //reset
						best_nets[net.Level()] = net;
						if(net.Level()>best_level_max) best_level_max=net.Level();
				}
				else {//failed to find a better net
					if(help.heuristics(net)>help.termination) failed_attempts[net.Level()]++;
					if (failed_attempts[net.Level()] % (help.max_fails/10) == 0 && failed_attempts[net.Level()]>0)
						cout << "! Failed attempt count at level "<<  net.Level() <<": " << failed_attempts[net.Level()] << "/" << help.max_fails <<"; open size=" << open.size() << endl;

					bool failed=true;
					for(int i=1;i<=best_level_max;i++) if (failed_attempts[i] < help.max_fails){ failed=false; break; }
					if(failed){
							cerr << "! Max level reached for all levels. Report the best unfolding." << endl;
							return best_nets[best_level_max];
					}
				}

				//determine how we should expand
				if( help.heuristics(net)>help.termination ) //expand when this node is popped from open for the first time
				{
					//expand to neighbors
					//cout<<"Expand in level="<<net.Level()<<endl;
					list<Net> neighbors = help.neighbors(net);

					//if(net.Level()==3) cout<<"nei size="<<neighbors.size()<<endl;

					for (Net& n : neighbors) {
						const BitVector& code = n.encode();
						//check if code has been visited
						//bool visited = (fscore.find(code)!=fscore.end());
						bool b = (visted.find(code) != visted.end());

						//since we know the distance from start, we don't update
						//the score (f value) ??
						if (!b) {
							//these two lines change the behave if deleted....
							//weird
							n.Level()=net.Level();
							help.dist(n,ecosts);
							help.heuristics(n);
							open.push_back(n);
							push_heap(open.begin(), open.end());
							visted.insert(code);
						}
					}//end for n
				}
				else{ //help.heuristics(net)<=help.termination

					if(net.Level()<help.max_level){
						//cout<<"Upgrade tp level="<<net.Level()+1<<endl;
						net.Level()++;
						net.G() = net.H() = FLT_MAX; //reset
						help.dist(net,ecosts); //recompute dist and heuristics
						help.heuristics(net);

						//cout<<"updated to level ="<<net.Level()<<" G=" << net.G()<<" H="<<net.H()<<endl;
						vector<Net>& open2=opens[current_level+1];
						open2.push_back(net);
						push_heap(open2.begin(), open2.end());
					}
					else{//solution found
						cout << "- Solution found! Open size=" << open.size() << " dist=" << help.dist(net,ecosts)
								 <<" h="<<help.heuristics(net)<< endl;
						return net;
					}
				}

				current_level=(current_level+1)%help.max_level;
			}//end while

			return best_nets[best_level_max];
		}
#endif

		//
		//
		// A* sub-modules
		//
		//

		list<EasilyFoldableUnfolding::Net>
		EasilyFoldableUnfolding::AStar_Helper_Foldable::neighbors(Net& net)
		{
			list<Net> neis;
			model * m=unfolder->getModel();
			int fbase=unfolder->getBaseFaceID();
			const BitVector& code=net.encode();

			//find edges connecting overlapping faces
			//these edges are to be cut!
			list<uint> edges_to_cut;

#ifdef LAZY
			auto m_cd_step_bkup=m_cd_step;
			m_cd_step/=pow(2, net.Level()-2);
			find_path(net, edges_to_cut);
			m_cd_step=m_cd_step_bkup;
#else
			find_path(net, edges_to_cut);
#endif

			//if(net.Level()==3) cout<<"edges_to_cut size="<<edges_to_cut.size()<<endl;

			//for (int i = 0; i < code.size(); i++) {
			for(uint i : edges_to_cut){
				const edge& e=m->edges[i];
				if (!code[i] || e.is2d() || e.type=='b') continue;
				list<int> nids = neighbor_edges2(code, m, i);
				for (int id : nids) {
					auto ncode = code;
					ncode.off(i);// = false;
					ncode.on(id);// = true;
					neis.push_back(Net(unfolder, fbase, ncode));
				}
			}//end for i

			return neis; //netghbors of this net
		}

		//-------------------
		//in this helper, we only flip the edges so there without affecting the keep faces

		float EasilyFoldableUnfolding::AStar_Helper_Foldable::dist
		(EasilyFoldableUnfolding::Net& net, const vector<float>& ecosts)
		{
				net.G()=0;
				return 0;

				model * m=unfolder->getModel();
				if (net.G() == FLT_MAX || g_recompute_dist) {
					float flips = 0;
					for (int i = 0; i < net.encode().size(); i++){
						if (this->root.encode()[i] != net.encode()[i]) {
							flips++;
							//const edge& e=m->edges[i];
							//flips+=e.length;
							//cout<<"Ahhhh..."<<endl;
						}
					}
					net.G()=flips;
				}

				return net.G();
		}

		float EasilyFoldableUnfolding::AStar_Helper_Foldable::heuristics
		(EasilyFoldableUnfolding::Net& net)
		{

				if (net.H() == FLT_MAX) {
					float h=0;
					//number of overlaps
					int n=unfolder->buildFromWeights(net.encode());

#ifdef LAZY
					if(net.Level()>1){
						auto m_cd_step_bkup=m_cd_step;
						m_cd_step/=pow(2, net.Level()-2);
						h=evaluate_folding_motion(unfolder);
						m_cd_step=m_cd_step_bkup;
					}
#else
					h=evaluate_folding_motion(unfolder);
#endif
					net.H()=h+n*1000;
				}

				return net.H();
		}


		//return the number of colliding face during the folding process
		int EasilyFoldableUnfolding::AStar_Helper_Foldable::evaluate_folding_motion(Unfolder * unfolder)
		{
					//cout<<"- Evaluate folding motion....."<<endl;
					unordered_set<int> colliding;
					model * mesh=unfolder->getModel();

					for(float i=1.0;i>=0;i-=m_cd_step)
					{
						unfolder->unfoldTo(i);
						if (i >= 0.990001){
							continue;
							/*
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
							*/
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

		//find a shortest path connecting f1 and f2
		//list<uint>& path contains a list of edges that the path passes through
		void EasilyFoldableUnfolding::AStar_Helper_Foldable::find_path(Net & net, list<uint>& path)
		{
			unfolder->buildFromWeights(net.encode());
			vector<set<uint> >  pairs;
			if(net.Level()==1)
				pairs=unfolder->getOverlppingFacePairs();
			else
					pairs=collect_self_intersection(net);

			//find a path connecting each pair of faces
			model * m=unfolder->getModel();
			set<uint> eid_set;
			for(int f1=0;f1<m->t_size;f1++){
					for(int f2 : pairs[f1]){
						list<uint> tmp;
						find_path(net, f1, f2, tmp );
						eid_set.insert(tmp.begin(), tmp.end());
					}
			}

			//the are set of folding edges that connect overlapping faces
			path.insert(path.end(), eid_set.begin(), eid_set.end());
		}

		void EasilyFoldableUnfolding::AStar_Helper_Foldable::find_path(Net & net, uint f1, uint f2, list<uint>& path)
		{

		  model * m=unfolder->getModel();

		  std::deque<uint> open;
		  bool * visited = new bool[m->t_size];
		  memset(visited, 0, sizeof(bool)* m->t_size);

		  //add root (f1)
		  visited[f1] = true;
		  open.push_back(f1);
		  m->tris[f1].parent_id = f1;

		  //loop until all faces in open are handled
			const BitVector& code=net.encode();
		  while (open.empty() == false)
		  {
		    uint fid = open.front();
		    open.pop_front();

		    if (fid == f2) //path found!
		    {
		      break; //done
		    }

		    for (int i=0;i<3;i++)
		    {
		      uint eid=m->tris[fid].e[i];
					edge & e=m->edges[eid];
		      //make sure eid is a crease edge
		      //if( net->is_crease(eid)==false ) continue;
					if( !code[eid] && e.is2d() ){
						cerr<<"! Error: edge "<<eid<<" is a 2d edge but is cut..."<<endl;
						exit(1);
					}
					if( !code[eid] ) continue;

		      uint kid=e.otherf(fid);
		      if (visited[kid]) continue;
		      open.push_back(kid);
		      visited[kid] = true;
		      m->tris[kid].parent_id = fid;
		    } //end for kid

		  } //end while

		  int now = f2;
		  while (now != f1)
		  {
		    int next=m->tris[now].parent_id;
		    path.push_front( m->getEdgeIdByFids(now,next) );
		    now = next;
		  }

		  delete [] visited;

		}


		//collect all intersections and place them inside the unfolder
		vector<set<uint>>
		EasilyFoldableUnfolding::AStar_Helper_Foldable::collect_self_intersection
		(const Net& net)
		{
				//return the number of colliding face during the folding process

				//cout<<"- Evaluate folding motion....."<<endl;
				unordered_set<int> colliding;
				model * mesh=unfolder->getModel();
				vector<set<uint>> collected_overlaps;
				collected_overlaps.resize(mesh->t_size);

				for(float i=1.0;i>=0;i-=m_cd_step)
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



		//get a list of edges if the eid is connected in n
		//list<int> BloomingUnfolding::Net::neighbor_edges(model * m, int eid)
		list<int> neighbor_edges2(const BitVector& code, model * m, int eid)
		{
			//compute connecected CCs after eid is cut
			DisjointSets dj(m->t_size);
			for (int i = 0; i < m->e_size; i++) {
				if (i == eid) continue; //ignore the edge to be cut
				edge& e = m->edges[i];
				if (e.type == 'b') continue;
				if (code[i] || e.is2d()) {
					int r1 = dj.find(e.fid.front());
					int r2 = dj.find(e.fid.back());
					if (r1 != r2) dj.unite(r1, r2);
				}
			}

			//get cut edges incident to these faces
			set<int> cut_edges;
			for (int i = 0; i < m->e_size; i++) {
					if (i == eid) continue; //ignore the edge to be cut
					edge& e = m->edges[i];
					if (e.is2d()) continue; //diagonal edge
					if (e.type == 'b') continue; //border edge
					if (code[i]) continue; //crease edge
					//now check if the faces connected by e are from different sets
					if (dj.find(e.fid.front()) != dj.find(e.fid.back()))
						cut_edges.insert(i);
			}//end for i

			///cout<<"cut edge size="<<cut_edges.size()<<endl;
			return list<int>(cut_edges.begin(), cut_edges.end());
		}

	}
}
