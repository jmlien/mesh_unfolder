/*
* BloomingUnfolding.cpp
*
*  Created on: April 22, 2020
*      Author: jmlien
*/
#include "BloomingUnfolding.h"

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

		//static data
		bool BloomingUnfolding::g_recompute_dist = false;
		short BloomingUnfolding::g_dist_level = 1;
		map<BloomingUnfolding::FACE_TYPE, float> BloomingUnfolding::g_face_cost =
		{ {BloomingUnfolding::BASE_FACE, 10000.0f},
		 {BloomingUnfolding::KEEP_FACE, 1.0f},
					 {BloomingUnfolding::TOSS_FACE, 0.001f} };

		//
		list<int> neighbor_edges(const BitVector& code, model * m, int eid);

		//get a list of  faces coplanar to the given face
		//the return list include the face itself
		//Note: this should actually be model to model.h
		inline void getCoplanarFaces(model * m, int face, list<int>& coplanars);

		BloomingUnfolding::BloomingUnfolding(Unfolder* unfolder)
		{
			const Config & mycfg = unfolder->getConfig();

			if (mycfg.find_best_base_face) {
				cerr << "! Error: BloomingUnfolding requires the flag -nbb" << endl;
				exit(1);
			}

			if (mycfg.baseface == -1) {
				cerr << "! Error: BloomingUnfolding requires the flag -bf to specify the base face" << endl;
				exit(1);
			}

			//this->m_unfolder = new Unfolder(unfolder->getModel(), mycfg);
			//assert(this->m_unfolder);

			//m_blooming_method = FLOWER;
			//m_splitter = NULL;
			//m_blooming_base_fase = mycfg.baseface;
			m_runs = mycfg.max_retries;
			m_start_time = 0;
			this->m_unfolder = unfolder;
			this->m_preview_only = false;
			this->m_no_trim=false;

			this->m_blooming_range=mycfg.blooming_range;
			if (mycfg.blooming_dir.normsqr() != 0) this->m_blooming_dir = mycfg.blooming_dir.normalize();
		}

		BloomingUnfolding::~BloomingUnfolding()
		{

		}

		bool BloomingUnfolding::setup(const string& method)
		{
/*
			if (method.find("flower") != string::npos) {
				cout << "- Blooming Unfolding uses flower ";
				m_blooming_method = FLOWER;
				m_splitter = new FlowerBloomingSplitter(m_unfolder->getModel(), m_blooming_base_fase);
			}
			else if (method.find("star") != string::npos) {
				cout << "- Blooming Unfolding uses star ";
				m_blooming_method = STAR;
				m_splitter = new StarBloomingSplitter(m_unfolder->getModel(), m_blooming_base_fase);
			}
			else {
				cerr << "! Error: BloomingUnfolding::setup: Unknow method: " << method << endl;
				return false;
			}
*/
			cout<<"(";
			if (method.find("preview") != string::npos)
			{
				this->m_preview_only = true; //enable preview only mode
				cout << "\"preview\" the initial unfolding only;";
			}
			else if (method.find("nocut") != string::npos)
			{
				this->m_no_trim = true; //avoid trimming the model
				cout << "\"preview\" A* unfolding only;";
			}
			else{
				cout << "add \"preview\" to see the initial unfolding;";
				cout << "add \"nocut\" to see the overlapping faces;";
			}
			cout<<")"<<endl;


			if(this->m_blooming_dir.normsqr()==0)
				this->m_blooming_dir = m_unfolder->getModel()->tris[m_unfolder->getBaseFaceID()].n;

			return true;
		}

		// will be called once after setup,
		// can be used to override parameters
		void BloomingUnfolding::init()
		{
			//nothing yet
		}

		void BloomingUnfolding::print(ostream& out) const
		{
			//nothing yet
		}


		bool BloomingUnfolding::blooming_unfold(Unfolder * unfolder)
		{
				auto mesh = unfolder->getModel();

				FlowerBloomingSplitter splitter(mesh, unfolder->getBaseFaceID());
				splitter.measure(mesh);

				Config mycfg = unfolder->getConfig();
				vector<float> best_weights;
				int min_overlap_count = INT_MAX;
				for (int r = 0; r < this->m_runs; r++)
				{
					//cfg.user_vector = rand_dirs[r];
					vector<float> weights = splitter.assignWeights(mesh, mycfg);
					auto count = unfolder->buildFromWeights(weights);
					if (count < min_overlap_count) {
						best_weights = weights;
						min_overlap_count = count;
						//cerr << "- iter = " << r << ", total overlaps = " << count << "\r" << flush;
						if (count == 0) break; //done, found a net
					}
				}
				//cerr<<"\n"<<flush;

				//finalize
				unfolder->buildFromWeights(best_weights);
				return min_overlap_count==0;
		}

		void BloomingUnfolding::run()
		{
			const auto config = this->m_unfolder->getConfig();
			model * mesh = this->m_unfolder->getModel();

			m_start_time = clock();

			//unfolding using blooming unfokding
			blooming_unfold(m_unfolder);


			//compute the face types from a given vector
			//in this caes the normal of the base face
			compute_face_types(m_unfolder);

			if (this->m_preview_only) return; //preview the unfolding without optimization

#define ASTAR 1
#if ASTAR == 3
			Net optmized_net = AStarMix();
#elif ASTAR == 2
			Net optmized_net = AStar();
#else
			Net optmized_net = AStar2Steps();
#endif

			m_unfolder->buildFromWeights(optmized_net.encode());

/*
			auto& overlaps = m_unfolder->getOverlppingFacePairs();
			for(int i=0;i<mesh->t_size;i++){
				if(overlaps[i].empty()) continue;
				cout<<"face["<<i<<"] intersects ";
				for(auto& o : overlaps[i]){
					cout<<o<<";";
				}
				cout<<endl;
			}
*/

			//remove faces that are overlapping....
			if(!this->m_no_trim){ //this is default
				cout<<"- Blooming unfolding finding optimal cuts to resolve overlappings"<<endl;
				const auto& overlaps=m_unfolder->getOverlppingFacePairs();
				optmized_net.optimalcuts(m_unfolder, overlaps);
				cout<<"- Blooming unfolding building new unfoldings"<<endl;
				m_unfolder = build_unfolder_from_net(m_unfolder, optmized_net);
				model * subm=m_unfolder->getModel();
				stringstream ss;
				ss << "./" << mesh->name << "_bloom_" << subm->t_size << ".obj";
				subm->saveObj(ss.str());
				cout<<"- Blooming unfolding saved the trimmed model to "<<ss.str()<<endl;
			}

			//int count=this->evaluate_folding_motion(m_unfolder);
			//cout<<"- collision count="<<count<<endl;

		}// END BloomingUnfolding::run()

		//same as AStar
		BloomingUnfolding::Net BloomingUnfolding::AStarMix()
		{
			//
			//get creases from the current unfolding
			//we assume that this is the best blooming unfolding
			//
			model * mesh=m_unfolder->getModel();
			const set<uint>& creases = m_unfolder->getFoldEdges();
			vector<bool> ic(mesh->e_size, 0);
			for (uint c : creases) ic[c] = true;

			//create a net from the creases
			Net start(m_unfolder, m_unfolder->getBaseFaceID(), BitVector(ic));
			AStar_Helper_Mixed mixed(start, m_unfolder, 1);

			return AStar(mixed);
		}

		//convert a current unfolding into a net
		BloomingUnfolding::Net BloomingUnfolding::getNet(Unfolder * unfolder)
		{
			model * m =unfolder->getModel();
			vector<bool> ic(m->e_size, 0);
			const set<uint>& creases = unfolder->getFoldEdges();
			for (uint c : creases) ic[c] = true;
			return Net(unfolder, unfolder->getBaseFaceID(), BitVector(ic));
		}

		//try to unfolde the keep faces and then toss faces
		BloomingUnfolding::Net BloomingUnfolding::AStar2Steps()
		{
			//this version create a submodel and unfold that submodel first
			model * mesh=m_unfolder->getModel();

			//rebuild unfolder by unfolding the KEEP_FACEs first
			//create a submodel with only KEEP faces that are connected to the base
			Unfolder * keep_unfolder = build_keep_unfolder(m_unfolder);
			//get creases from the current unfoldings
			model * keep_mesh=keep_unfolder->getModel();
			keep_mesh->saveObj("keep_mesh.obj");

			{//find optimal unfolding for this keep_unfolder
				Net start=getNet(keep_unfolder);
				compute_face_types(keep_unfolder);
				AStar_Helper_Mixed keep(start, keep_unfolder, 0,3,1000); //stop when there is 0 overlap and try it until level_3
				Net keep_net = AStar(keep);
				keep_unfolder->buildFromWeights(keep_net.encode());

				//this->m_unfolder=keep_unfolder;
				//return keep_net;
			}


			{//rebuild m_unfoder using keep_unfolder
				vector<float> ic(mesh->e_size, FLT_MAX);
				{//get data from the first unfolder
					const set<uint>& creases = m_unfolder->getFoldEdges();
					for (uint c : creases) ic[c] = 1;
				}
				{//transfer data from the keep_unfolder
					const set<uint>& creases = keep_unfolder->getFoldEdges();
					for (uint c : creases){
						const edge& e=keep_mesh->edges[c];
						assert(e.type!='d');
						const triangle& t1=keep_mesh->tris[e.fid.front()];
						const triangle& t2=keep_mesh->tris[e.fid.back()];
						int eid=mesh->getEdgeIdByFids(t1.source_fid,t2.source_fid);
						assert(eid>=0 &&eid<mesh->e_size);
						ic[eid] = 0.1;
					}
				}

				//rebuild from this weight
				m_unfolder->buildFromWeights(ic);
			}//done rebuild keep_unfolder

			//
			Net start=getNet(m_unfolder);
			start.compute_face_types(m_unfolder); //update the face type so everything that is not
			AStar_Helper_Toss toss(start, m_unfolder, 0, 1,1000); //stop when there is 0 overlap and try it only for level_1
			return AStar(toss);
		}

		BloomingUnfolding::Net BloomingUnfolding::AStar(AStar_Helper& help)
		{
			//build the first net from unfolder
			Unfolder * unfolder=help.unfolder;
			model * m = unfolder->getModel();
			BloomingUnfolding::g_dist_level = 1; //start with level 1: full distance
			const int fbase = unfolder->getBaseFaceID(); //base face of a net, same of all nets

			//compute additional data to help guiding A* search
			vector<float> ecosts; //cost of flipping a root edge
			help.root.compute_edge_cost(m, ecosts);

			//start A* seach from the start net
			vector<Net> open;
			set<BitVector> visted;

			Net best_net = help.root;
			visted.insert(help.root.encode());
			open.push_back(help.root);

			int failed_attempts = 0;

			cout << "- Start A* optimization"<<endl;
			while (!open.empty()) {

				Net net = open.front();
				pop_heap(open.begin(), open.end());
				open.pop_back();

				///updating the conditions
				if (help.heuristics(net)<=help.termination) //research the goal
				{
					cout << "- Solution found! Open size=" << open.size() << " dist=" << help.dist(net,ecosts)
					     <<" h="<<help.heuristics(net)<< endl;
					return net;
				}
				else if ( (help.heuristics(net)+1 <= help.heuristics(best_net)) ||
				( (help.heuristics(net)== help.heuristics(best_net)) && help.dist(net,ecosts)<help.dist(best_net,ecosts)) )
				{
						cout << "- Best net overlaps=" << help.heuristics(net)
							<< " dist=" << help.dist(net,ecosts) << " open size=" << open.size() << endl;
						//reset only when the overall score is descreased
						if (help.heuristics(net)+help.dist(net,ecosts) <= help.heuristics(best_net)+help.dist(best_net,ecosts))
							failed_attempts = 0; //reset
						best_net = net;
				}
				else {//failed to find
					failed_attempts++;
					if (failed_attempts % (help.max_fails/10) == 0)
						cout << "! Failed attempt count: " << failed_attempts << "/" << help.max_fails << endl;

					if (failed_attempts > help.max_fails) {

						if (BloomingUnfolding::g_dist_level == help.max_level) { //still failed as max level?
							cerr << "! Max level reached. Report the best unfolding." << endl;
							return best_net;
						}

						cout << "- Max failed attempts reached (" << help.max_fails << "); change distance metrics to level " << BloomingUnfolding::g_dist_level + 1
							<< " open size=" << open.size() << endl;

						//OK, move to the next level of distance metric
						BloomingUnfolding::g_dist_level++;
						BloomingUnfolding::g_recompute_dist = true;
						for (Net& net : open) help.dist(net,ecosts); //update the distance
						//help.dist(best_net,ecosts);
						BloomingUnfolding::g_recompute_dist = false;
						make_heap(open.begin(), open.end()); //since distance changed

						failed_attempts = 0; //reset the counter
					}
				}

				//expand to neighbors
				list<Net> neighbors = help.neighbors(net);
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
						float g = help.dist(n,ecosts);
						float f = g + help.heuristics(n);
						open.push_back(n);
						push_heap(open.begin(), open.end());
						visted.insert(code);
					}
				}//end for n

			}//end while

			return best_net;
		}






		BloomingUnfolding::Net BloomingUnfolding::AStar()
		{
			//build the first net from unfolder
			model * m = m_unfolder->getModel();
			BloomingUnfolding::g_dist_level = 1; //start with level 1: full distance

			//get creases from the current unfolding
			//we assume that this is the best blooming unfolding
			const set<uint>& creases = m_unfolder->getFoldEdges();
			vector<bool> ic(m->e_size, 0);
			for (uint c : creases) ic[c] = true;

			//create a net from the creases
			const int fbase = m_unfolder->getBaseFaceID(); //base face of a net, same of all nets
			Net start(m_unfolder, fbase, BitVector(ic));

			//compute addition data to help guiding A* search
			vector<float> ecosts; //cost of flipping a root edge
			compute_face_types(m_unfolder);
			start.compute_edge_cost(m, ecosts);

			//start A* seach from the start net
			vector<Net> open;
			Net best_net = start;

			//map<BitVector, float> gscore;
			//map<BitVector, float> fscore;
			//set<BitVector> visted;
			unordered_set<BitVector,hash_BitVector> visted;

			//set<BitVector> closed;
			//unordered_set<BitVector,hash_BitVector> closed;
			//fscore[start.encode()] = start.heuristics(m_unfolder);
			visted.insert(start.encode());
			open.push_back(start);


			const int max_failed_attempts = 1000;
			int failed_attempts = 0;

			cout << "- Start A* optimization with max failed attempt = " << max_failed_attempts << endl;
			while (!open.empty()) {
				Net net = open.front();
				pop_heap(open.begin(), open.end());
				open.pop_back();

				if (net.heuristics(m_unfolder) < 1) //research the goal
				{
					cout << "- Solution found! Open size=" << open.size() << " dist=" << net.dist(m,start,ecosts)
					     <<" h="<<net.heuristics(m_unfolder)<< endl;
					/*
										vector<bool> value;
										net.encode().tovector(value);
										net.analyze(m, fbase, value);
					*/
					return net;
				}
				else if (net.heuristics(m_unfolder) < best_net.heuristics(m_unfolder))
				{
					//if(net.heuristics(m_unfolder) >= 1)
					{
						best_net = net;
						cout << "- Best net overlaps=" << net.heuristics(m_unfolder)
							<< " dist=" << net.dist(m,start,ecosts) << " open size=" << open.size() << endl;
						failed_attempts = 0; //reset
					}
				}
				else {//failed to find
					failed_attempts++;
					if (failed_attempts % (max_failed_attempts/10) == 0)
						cout << "! Failed attempt count: " << failed_attempts << "/" << max_failed_attempts << endl;

					// if(net.heuristics(m_unfolder) < 1) //only non-essential faces overlap
					// {
					// 	//already try many times or already in level3
					// 	if(failed_attempts > max_failed_attempts || BloomingUnfolding::g_dist_level == 2)
					// 	{
					// 		cerr<<"! Max level reached. Only non-essential faces overlap. Move on."<<endl;
					// 		return best_net;
					// 	}
					// }

					if (failed_attempts > max_failed_attempts) {

						if (BloomingUnfolding::g_dist_level == 3) { //still failed as max level?
							cerr << "! Max level reached. Report the best unfolding." << endl;
							return best_net;
						}

						cout << "- Max failed attempts reached (" << max_failed_attempts << "); change distance metrics to level " << BloomingUnfolding::g_dist_level + 1
							<< " open size=" << open.size() << endl;

						//OK, move to the next level of distance metric
						BloomingUnfolding::g_dist_level++;
						BloomingUnfolding::g_recompute_dist = true;
						for (Net& net : open) net.dist(m,start,ecosts); //update the distance
						BloomingUnfolding::g_recompute_dist = false;
						make_heap(open.begin(), open.end()); //since distance changed
						// if(open.size()>max_open_size){
						// 	vector<Net> tmp;
						// 	for(int i=0;i<max_open_size/2;i++){
						// 		tmp.push_back(open.front());
						// 		open.pop
						// 	}
						// }
						failed_attempts = 0; //reset the counter
					}
				}

				//closed.insert(net.encode());

				//expand to neighbors
				list<Net> neighbors = net.neighbors(m_unfolder, fbase);
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
						float g = n.dist(m,start,ecosts);
						float f = g + n.heuristics(m_unfolder);
						open.push_back(n);
						push_heap(open.begin(), open.end());

						//if(open.size()>max_open_size) open.pop_back();

						//fscore[code]=f;
						visted.insert(code);
					}
				}//end for n

			}//end while

			return best_net;
		}

		//build a new unfolder from the net, this is used most likely that
		//some faces should be removed to ensure that there is no overlapping
		Unfolder * BloomingUnfolding::build_unfolder_from_net
		(Unfolder *unfolder, BloomingUnfolding::Net& net)
		{
			Config cfg = unfolder->getConfig();
			BitVector B = net.encode();
			model * m = unfolder->getModel();
			model * subm = build_model_from_net(m, unfolder->getBaseFaceID(), B);

			cfg.baseface = 0; //now the base face is 0
			Unfolder * sub_unfolder = (subm == m) ? unfolder : new Unfolder(subm, cfg);
			// vector<bool> values;
			// B.tovector(values);
			sub_unfolder->buildFromWeights(B);
			return sub_unfolder;
		}


		//build an unfolder that has only the KEEP_FACEs
		Unfolder * BloomingUnfolding::build_keep_unfolder(Unfolder *unfolder)
		{
			Net net=getNet(unfolder);
			model * m = unfolder->getModel();
			Config cfg = unfolder->getConfig();

			//collect faces that are connected to the root
			BitVector B(vector<bool>(m->e_size,true));
			int fbase=unfolder->getBaseFaceID();
			list<uint> fids = getCC(m,fbase,B,TOSS_FACE);

			//build a submodel from this list of faces
			model *subm=build_submodel(m,fids);

			//find the new base id...
			cfg.baseface=0;
			Unfolder * sub_unfolder = new Unfolder(subm, cfg);

			//unfold
			blooming_unfold(sub_unfolder);

			//done
			return sub_unfolder;
		}

		//build a model from the bit vector
		//this likely is used when some faces of the input net is trimmed by the ned
		//the bitvector may be updated to reflect the new model
		model * BloomingUnfolding::build_model_from_net(model * m, int base, BitVector& B)
		{
			//we will need a list of faces first, if the number if the same as m->t_size
			//then we are done.
			list<uint> fids = getCC(m,base,B);
			model * subm=build_submodel(m, fids);

			//rebuild B using the new ids
			BitVector B2(subm->e_size);
			for (int eid = 0; eid < subm->e_size; eid++)
			{
				edge& e = subm->edges[eid];
				if (e.type == 'b') continue;

				int fid1 = subm->tris[e.fid.front()].source_fid;
				int fid2 = subm->tris[e.fid.back()].source_fid;
				int old_id = m->getEdgeIdByFids(fid1, fid2);
				if (old_id < 0) {
					cerr << "! Error: BloomingUnfolding::build_submodel_from_net: "
						<< "Cannot find e" << endl;
					exit(1);
				}
				if (B[old_id]) B2.on(eid);
			}
			B = B2;
			return subm;
		}


		//create a subset of m using the fids
		model * BloomingUnfolding::build_submodel(model * m, list<uint> fids)
		{
			masc::obj::objModel data;

			//collect all vertices
			unordered_map<uint, uint> vids;
			for (auto fid : fids)
			{
				for (short d = 0; d < 3; d++) {
					vids[m->tris[fid].v[d]] = 0;
				}
			}

			//get a list of points
			uint new_vid = 0;
			for (auto & vid : vids)
			{
				masc::obj::Vpt pt;
				auto & pos = m->vertices[vid.first].p;
				pt.x = pos[0];
				pt.y = pos[1];
				pt.z = pos[2];
				data.pts.push_back(pt);
				vid.second = new_vid++;
			}

			//get a list of faces
			for (auto fid : fids)
			{
				masc::obj::polygon poly;
				for (short d = 0; d < 3; d++) {
					auto vid = vids[m->tris[fid].v[d]];
					poly.pts.push_back(vid);
				}
				data.polys.push_back(poly);
			}

			data.compute_v_normal();

			//build mesh
			model * subm = new model();
			subm->name=m->name;
			if (subm->build(data, true) == false)
			{
				cerr << "! Error: Failed to build a model from net" << endl;
				return NULL;
			}

			subm->makeManifold();

			//register the new fid to the old fid
			auto fid_it = fids.begin();
			for (int fid = 0; fid < subm->t_size; fid++)
			{
				subm->tris[fid].source_fid = *fid_it;
				++fid_it;
			}

			return subm;
		}

		void BloomingUnfolding::compute_face_types(Unfolder * unfolder)
		{
			compute_face_types(unfolder, m_blooming_dir, m_blooming_range);
		}

		//compute the type of each face in the model using the given viewing dirsection
		//store face type in data "cluster_id"
		void BloomingUnfolding::compute_face_types
		(Unfolder * unfolder, const Vector3d& dir, float range)
		{
			model * m = unfolder->getModel();

			for (int i = 0; i < m->t_size; i++) {
				triangle& tri = m->tris[i];
				if (i == unfolder->getBaseFaceID()) { tri.cluster_id = BASE_FACE; }
				else {
					//we abuse cluster_id here
					tri.cluster_id = (dir*tri.n >= -range) ? KEEP_FACE : TOSS_FACE;
				}//end if

				//cout<<"t["<<  i<<"] type="<<tri.cluster_id<<endl;
			}//end for i
		}


				list<uint> BloomingUnfolding::getCC(model * m, int base, const BitVector& B, int toss_face_type)
				{
					//we will need a list of faces first, if the number if the same as m->t_size
					//then we are done.
					list<uint> fids;
					{//this should be moved to a function actually
						list<int> open;
						open.push_back(base);
						vector<bool> visited(m->t_size, false);
						visited[base] = true;

						while (!open.empty()) {
							int n = open.front();
							open.pop_front();
							fids.push_back(n); //remember all the faces we visited


							for (int i = 0; i < 3; i++) {
								int eid = m->tris[n].e[i];

								const edge& e = m->edges[eid];
								if (!B[eid] && e.type != 'd') continue; //not connected
								if (e.type == 'b') continue; //there is no otherf
								int of = e.otherf(n);
								if (visited[of]) continue; //visited
								visited[of]=true;
								const triangle& t=m->tris[of];
								if( t.cluster_id == toss_face_type ) continue; //toss away
								open.push_back(of);
							}//edn for
						}//end while
					}
					//return list<uint>(fids.begin(),fids.end());
					return fids;
				}


//
//
// A* sub-modules
//
//

//-------------------
//mixed all faces and optimize them together
list<BloomingUnfolding::Net>
BloomingUnfolding::AStar_Helper_Mixed::neighbors
(BloomingUnfolding::Net& net)
{
	list<Net> neis;
	model * m=unfolder->getModel();
	int fbase=unfolder->getBaseFaceID();
	const BitVector& code=net.encode();

	for (int i = 0; i < code.size(); i++) {
		const edge& e=m->edges[i];
		if (!code[i] || e.type=='d' || e.type=='b') continue;
		list<int> nids = neighbor_edges(code, m, i);
		for (int id : nids) {
			auto ncode = code;
			ncode.off(i);// = false;
			ncode.on(id);// = true;
			neis.push_back(Net(unfolder, fbase, ncode));
		}
	}//end for i

	return neis; //netghbors of this net
}

float BloomingUnfolding::AStar_Helper_Mixed::dist
(BloomingUnfolding::Net& net, const vector<float>& ecosts)
{
	if (net.G() == FLT_MAX || g_recompute_dist) {

		float g=0;
		model * m=unfolder->getModel();

		if (g_dist_level <= 2) {

			//# of bits flipped
			float flips = 0;
			for (int i = 0; i < m->e_size; i++){
				if (this->root.encode()[i] != net.encode()[i]){
					const edge& e=m->edges[i];
	#if 0
					//flips++;
					flips+=e.length;
					//flips+=ecosts[i]/10000.0f;
	#else
					const triangle& t1=m->tris[e.fid.front()];
					const triangle& t2=m->tris[e.fid.back()];
					auto type_1=(BloomingUnfolding::FACE_TYPE)t1.cluster_id;
					auto cost_1=BloomingUnfolding::g_face_cost.at(type_1);
					auto type_2=(BloomingUnfolding::FACE_TYPE)t2.cluster_id;
					auto cost_2=BloomingUnfolding::g_face_cost.at(type_2);
					flips+=(cost_1*cost_2); //e.length;
	#endif
				}
			}

			//float area_score = (net.getHullArea()-root.getHullArea())/root.getHullArea();
			g = flips; // + area_score;

			//done
			if (g_dist_level == 1) {

				//decrease of # of leaves
				int leaf_score = (root.getLeafSize() - net.getLeafSize()); //more leaf is better
				//if (leaf_score < 0) leaf_score = 0;

				//increase of diameter
				int dia_score = (net.getDepth() - root.getDepth()); //weighted by 100
				//if (dia_score < 0) dia_score = 0;

				//cout<<"area_score="<<area_score<<endl;

				g += leaf_score + dia_score;
			}
		}
		else g = 0; //g_dist_level>2

		net.G()=g;
	}

	return net.G();
}

float BloomingUnfolding::AStar_Helper_Mixed::heuristics
(BloomingUnfolding::Net& net)
{

				if (net.H() == FLT_MAX) {
					float h=0;
					//number of overlaps
					// vector<bool> creases;
					// this->code.tovector(creases);
					int n=unfolder->buildFromWeights(net.encode());
	#if 0
					const vector<set<uint>>& overlaps = unfolder->getOverlppingFacePairs();
					int overlap_count = 0;
					for(auto& o : overlaps ) if(!o.empty()) overlap_count++;
					h = overlap_count;
					//h=n;
	#else
					h = overlapping_cost(unfolder);
	#endif

					net.H()=h;
				}

				return net.H();
}

//-------------------
//in this helper, we only flip the edges without affecting the keep faces
list<BloomingUnfolding::Net> BloomingUnfolding::AStar_Helper_Toss::neighbors(BloomingUnfolding::Net& net)
{
	list<Net> neis;
	model * m=unfolder->getModel();
	int fbase=unfolder->getBaseFaceID();
	const BitVector& code=net.encode();

	//go through each edge
	for (int i = 0; i < m->e_size; i++) {
		const edge& e=m->edges[i];
		if (!code[i] || e.type=='d' || e.type=='b') continue;
		if(isKeep(i)) continue; //cannot be a keep edge

		//if (Keep_SUBNET(net, i)) continue;

		list<int> nids = neighbor_edges(code, m, i);
		for (int id : nids) {
			if(isKeep(id)) continue; //cannot use keep edge

			//if (Keep_SUBNET(net, id)) continue;

			auto ncode = code;
			ncode.off(i);// = false;
			ncode.on(id);// = true;
			neis.push_back(Net(unfolder, fbase, ncode));
		}
	}//end for i

	return neis; //netghbors of this net
}

float BloomingUnfolding::AStar_Helper_Toss::dist(BloomingUnfolding::Net& net, const vector<float>& ecosts)
{
	if (net.G() == FLT_MAX || g_recompute_dist) {

		float g=0;
		model * m=unfolder->getModel();

		if (g_dist_level <= 2) {

			//# of bits flipped
			float flips = 0;
			for (int i = 0; i < m->e_size; i++){
				if (this->root.encode()[i] != net.encode()[i]){
					if(isKeep(i)){
						cerr<<"! Error: BloomingUnfolding::AStar_Helper_Toss::dist: Edges should NOT be KEEP edges"<<endl;
						exit(1);
					}
					flips++;
				}
			}
			g = flips;

			//done
			if (g_dist_level == 1) {
				//decrease of # of leaves
				int leaf_score = (root.getLeafSize() - net.getLeafSize()); //more leaf is better
				//if (leaf_score < 0) leaf_score = 0;

				//increase of diameter
				int dia_score = (net.getDepth() - root.getDepth()); //weighted by 100
				//if (dia_score < 0) dia_score = 0;

				g += leaf_score + dia_score;
			}
		}
		else g = 0; //g_dist_level>2

		net.G()=g;
	}

	return net.G();
}

float BloomingUnfolding::AStar_Helper_Toss::heuristics(BloomingUnfolding::Net& net)
{
	if (net.H() == FLT_MAX) {
		float h=0;

		//unfold
		model * m=unfolder->getModel();
		unfolder->buildFromWeights(net.encode());

		const vector<set<uint>>& overlaps = unfolder->getOverlppingFacePairs();
		int tossxtoss_count = 0;
		int tossxkeep_count = 0;
		int keepxkeep_count = 0;
		for(int i=0;i<m->t_size;i++){
			const set<uint>& o=overlaps[i];
			if(o.empty()) continue;
			const triangle& ti=m->tris[i];
			auto type_i=(BloomingUnfolding::FACE_TYPE)ti.cluster_id;

			if (type_i == BloomingUnfolding::BASE_FACE) {//VERY BAD, no need to look further
				return net.H() = FLT_MAX / 2;
			}

			for(int j : o){
				if (i > j) continue;
				const triangle& tj=m->tris[j];
				auto type_j=(BloomingUnfolding::FACE_TYPE)tj.cluster_id;
				if(type_i==type_j && type_i==BloomingUnfolding::KEEP_FACE){
					//cerr<<"! Warning: BloomingUnfolding::AStar_Helper_Toss::heuristics: There should not be any KEEP face overlap"<<endl;
					keepxkeep_count++;
				}
				if(type_i==type_j && type_i==BloomingUnfolding::TOSS_FACE)
					tossxtoss_count++;
				else
					tossxkeep_count++;
			}//end for j
		}//end for i

	#if 0
		{
			h = tossxkeep_count+tossxtoss_count;
		}
	#else
		{
			h = keepxkeep_count*100+tossxkeep_count+tossxtoss_count*0.01;
			/*
			auto toss_bkup=g_face_cost[TOSS_FACE];
			auto keep_bkup=g_face_cost[KEEP_FACE];
			g_face_cost[TOSS_FACE]=0.001; //
			g_face_cost[KEEP_FACE]=1000;  //this makes keepxkeep 1e6 and keepxtoss 1, so we can have at most 1e6-1 pairs of keepxtoss overlaps
			h = overlapping_cost(unfolder);
			g_face_cost[TOSS_FACE]=toss_bkup;
			g_face_cost[KEEP_FACE]=keep_bkup;
			*/
		}
	#endif

		net.H()=h;
	}

	return net.H();
}

bool BloomingUnfolding::AStar_Helper_Toss::isKeep(int eid) //is the edge incident to keep faces
{
	model * m=unfolder->getModel();
	const edge& e=m->edges[eid];

	for(int f : e.fid){
		const triangle& t = m->tris[f];
		auto type=(BloomingUnfolding::FACE_TYPE)t.cluster_id;
		if(type!=BloomingUnfolding::KEEP_FACE) return false;
	}

	return true;
}


		  //
			//
			// Net methods
			//
			//

		BloomingUnfolding::Net::Net(const BloomingUnfolding::Net & other) //copy constructor
			:code(other.code)
		{
			this->root = other.root;
			//this->creases = other.creases;
			this->g = other.g;
			this->h = other.h;
			this->n = other.n;
			this->level=other.level;
			this->depth = other.depth;
			this->leaves = other.leaves;
			this->parent = other.parent;
			this->kids = other.kids;
			this->m_overlaps = other.m_overlaps; //overlapping pairs
			this->m_hull_area=other.m_hull_area ; //convex hull area
		}

		// BloomingUnfolding::Net::Net(model * m, int root, const vector<bool>& value) //number of edges in this net
		// 	:BloomingUnfolding::Net::Net(m,root,BitVector(value))
		// {}

		BloomingUnfolding::Net::Net
		(Unfolder * unfolder, int root, const BitVector& value) //number of edges in this net
			:code(value)
		{
			//init data
			this->root = root;
			this->level=1;
			g = h = m_hull_area=FLT_MAX;
			this->n = code.size();
			analyze(unfolder, code);
		}

		//get the neighbors of this net
		list<BloomingUnfolding::Net> BloomingUnfolding::Net::neighbors
		(Unfolder * unfolder, int fbase)
		{
			list<Net> neis;
			model * m=unfolder->getModel();
			for (int i = 0; i < n; i++) {
				const edge& e=m->edges[i];
				if (!this->code[i] || e.type=='d' || e.type=='b') continue;
				list<int> nids = neighbor_edges(this->code, m, i);
				for (int id : nids) {
					auto ncode = this->code;
					ncode.off(i);// = false;
					ncode.on(id);// = true;

					// {
					  //   cout << "-----------------" << endl;
					  //   int sum = 0;
					  //   for (bool b : ncreases) sum += b;
					  //   if (sum != m->t_size - 1) {
					  // 	  cout << "number of connected edge is  " << sum << " but should be " << m->t_size - 1 << endl;
					  // 	  exit(1);
					  //   }
					// }

					neis.push_back(Net(unfolder, fbase, ncode));
				}
			}//end for i

			return neis; //netghbors of this net
		}

		//the distance from the root
		float BloomingUnfolding::Net::dist(model * m, const Net& root, vector<float>& ecosts)
		{
			if (g == FLT_MAX || g_recompute_dist) {

				if (g_dist_level <= 2) {
					//# of bits flipped
					float flips = 0;
					for (int i = 0; i < this->n; i++){
						if (root.code[i] != this->code[i]){
							const edge& e=m->edges[i];
#if 1
							flips++;
							//flips+=ecosts[i];
#else
							const triangle& t1=m->tris[e.fid.front()];
							const triangle& t2=m->tris[e.fid.back()];
							auto type_1=(BloomingUnfolding::FACE_TYPE)t1.cluster_id;
							auto cost_1=BloomingUnfolding::g_face_cost.at(type_1);
							auto type_2=(BloomingUnfolding::FACE_TYPE)t2.cluster_id;
							auto cost_2=BloomingUnfolding::g_face_cost.at(type_2);
							flips+=(cost_1*cost_2); //e.length;
#endif
						}
					}
					g = flips;

					//done
					if (g_dist_level == 1) {
						//decrease of # of leaves
						int leaf_score = (root.leaves - this->leaves); //more leaf is better
						//if (leaf_score < 0) leaf_score = 0;


						//increase of diameter
						int dia_score = (this->depth - root.depth); //weighted by 100
						//if (dia_score < 0) dia_score = 0;

						g += leaf_score + dia_score;
					}
				}
				else g = 0; //g_dist_level>2
			}

			return g;
		}

		//distance to go heuristics
		float BloomingUnfolding::Net::heuristics(Unfolder* unfolder)
		{

			if (h == FLT_MAX) {
				//number of overlaps
				// vector<bool> creases;
				// this->code.tovector(creases);
				int n=unfolder->buildFromWeights(this->code);
#if 0
				const vector<set<uint>>& overlaps = unfolder->getOverlppingFacePairs();
				int overlap_count = 0;
				for(auto& o : overlaps ) if(!o.empty()) overlap_count++;
				h = overlap_count;
				//h=n;
#else
				h = overlapping_cost(unfolder);
#endif
			}

			return h;
		}

		//get the hash code of this net
		const BitVector& BloomingUnfolding::Net::encode() const
		{
			return this->code;
		}

		//comparator
		bool BloomingUnfolding::Net::operator<(const BloomingUnfolding::Net & other) const
		{
			if (h == FLT_MAX || g == FLT_MAX || other.h== FLT_MAX || other.g== FLT_MAX){
					cerr<<"! Error: BloomingUnfolding::Net::operator<: Net is not initialized"<<endl;
					cout<<"h="<<h<<" g="<<g<<" other.h="<<other.h<<" other.g="<<other.g<<endl;
					exit(1);
			}

			float f = this->h + this->g;
			float of = other.h + other.g;

			if (f != of)  return f > of; //">" makes our heap a min heap
			return this->h > other.h;

			//if(this->h != other.h) return this->h > other.h;
			//else{
			//  return this->h + this->g>other.h + other.g;
			//}
		}

		//get a list of edges if the eid is connected in n
		//list<int> BloomingUnfolding::Net::neighbor_edges(model * m, int eid)
		list<int> neighbor_edges(const BitVector& code, model * m, int eid)
		{
			//compute connecected CCs after eid is cut
			DisjointSets dj(m->t_size);
			for (int i = 0; i < m->e_size; i++) {
				if (i == eid) continue; //ignore the edge to be cut
				edge& e = m->edges[i];
				if (e.type == 'b') continue;
				if (code[i] || e.type == 'd') {
					int r1 = dj.find(e.fid.front());
					int r2 = dj.find(e.fid.back());
					if (r1 != r2) dj.unite(r1, r2);
				}
			}

			//get faces from adjacent faces
			list<int> adj_faces;
			for (int fid : m->edges[eid].fid)
				getCoplanarFaces(m, fid, adj_faces);

			//get cut edges incident to these faces
			//JML: Why? Can't we use all cut edges?
			set<int> cut_edges;
			for (int fid : adj_faces) {
				triangle& tri = m->tris[fid];
				for (int i = 0; i < 3; i++) {
					int x = tri.e[i];
					if (x == eid) continue;
					edge& e = m->edges[x];
					if (e.type == 'd') continue; //diagonal edge
					if (e.type == 'b') continue; //border edge
					if (code[x]) continue; //crease edge

					//now check if the faces connected by e are from different sets
					if (dj.find(e.fid.front()) != dj.find(e.fid.back()))
						cut_edges.insert(x);
				}//end for i
			}//end for fid

			return list<int>(cut_edges.begin(), cut_edges.end());
		}


		//compute some statistics of this net, directly from code
		//0. build the tree
		//1. count # of leaves
		//2. distance from the base to the farthest leaf
		void BloomingUnfolding::Net::analyze
		(Unfolder * unfolder, const BitVector& value)
		{
			model * m =unfolder->getModel();
/*
			auto count = unfolder->buildFromWeights(value);

			//collecte intersecting pairs....
			const vector<set<uint>>& overlaps = unfolder->getOverlppingFacePairs();
			for(uint i=0;i<m->t_size;i++)
				for(uint j: overlaps[i])
					if(i<j) this->m_overlaps.push_back(make_pair(i,j));
*/
			//compute the areas
			// this->m_hull_area=unfolder->getHullArea();

			//cout<<"this->m_hull_area="<<this->m_hull_area<<endl;

			//---------------------------------
			list<int> open;
			open.push_back(root);
			vector<int> depth(m->t_size, -1);
			depth[root] = 0;
			int leaf_count = 0;
			int max_depth = 0;
			this->parent.resize(m->t_size);      //parent[i] is the parent of face i
			this->kids.resize(m->t_size);        //kids[i] is the kids of face i
			parent[root] = root;

			while (!open.empty()) {
				int n = open.front();
				open.pop_front();
				if (max_depth < depth[n]) max_depth = depth[n];

				int kid_count = 0; //# of kids
				for (int i = 0; i < 3; i++) {
					int eid = m->tris[n].e[i];

					edge& e = m->edges[eid];
					if (!value[eid] && e.type != 'd') continue; //not connected

					int of = e.otherf(n);
					if (depth[of] >= 0) continue; //visited

					depth[of] = depth[n] + 1;
					open.push_back(of);
					kids[n].push_back(of);
					parent[of] = n;
					kid_count++;
				}
				if (kid_count == 0) leaf_count++; //no kid, this is a leaf
			}

			//debug
				  //for (int d : depth) if (d == -1) { cerr << "ERROR depth; net has more than one cc" << endl; exit(1); }

			this->leaves = leaf_count;
			this->depth = max_depth;
		}

		//rebuild this net from the give code
		void BloomingUnfolding::Net::rebuild(Unfolder * unfolder, const BitVector& B)
		{
			//init data
			this->code = B;
			this->depth = 0;
			this->leaves = 0;
			this->parent.clear();
			this->kids.clear();
			this->m_overlaps.clear();

			m_hull_area=FLT_MAX;
			g = h = FLT_MAX;
			this->n = (int)this->code.size();
			analyze(unfolder, this->code);
		}

		//naive method
		//compute a list of areas that is in the subnet of a give face
		//this value include the area of the face itself
		//pre consition: vector<float>& areas must contain m->t_size zeros when this function is first
		//called
		float BloomingUnfolding::overlapping_cost(Unfolder *unfolder)
		{
			const vector<set<uint>>& overlaps = unfolder->getOverlppingFacePairs();
			model * m = unfolder->getModel();
			float cost = 0;

			//multiplitive
			for (uint i = 0; i < m->t_size; i++) {//for each face
				if (overlaps[i].empty()) continue; //nothing here, move on
				const triangle& ti = m->tris[i];
				if (ti.cluster_id == BloomingUnfolding::BASE_FACE) return FLT_MAX / 2;
				auto type_i=(BloomingUnfolding::FACE_TYPE)ti.cluster_id;
				auto cost_i=BloomingUnfolding::g_face_cost.at(type_i);
				for (uint j : overlaps[i]) {
					if (i > j) continue;
					const triangle& tj = m->tris[j];
					if (tj.cluster_id == BloomingUnfolding::BASE_FACE) return FLT_MAX / 2;
					//i overlaps with j
					auto type_j=(BloomingUnfolding::FACE_TYPE)tj.cluster_id;
					auto cost_j=BloomingUnfolding::g_face_cost.at(type_j);
					cost += cost_i*cost_j;
					//cost++;
				}
			}//end for i

			return cost;
		}


		//find optimal trim (instead of cut) of the net to avoid collision
		//in this case, we will remove the overlapping part of the face instead of cutting off the entire face
		//the code that trime the face will be performed later in method: XXXX
		//the input "overlap" would be modified to record the faces to be trimmed
		//note that in this case, the net won't change; it is only the mesh that will need to be modified
		float BloomingUnfolding::Net::optimaltrims(Unfolder *unfolder, vector<set<uint>>& overlaps)
		{
			model * m = unfolder->getModel();
			BitVector killed(m->t_size); //record what face is removed
			BitVector new_code=this->code;

			vector<float> face_areas(m->t_size, 0); //
			compute_area_in_subnet(m, this->root, face_areas);
			list<int> Q; //a queue
			Q.push_back(unfolder->getBaseFaceID());
			float cost = optimalcuts(unfolder, Q, new_code, killed, face_areas, overlaps);

			//TODO:
			//compute the difference between the new_code and this->code

			//the difference will tell us where the optimal cuts would be made
			//update overlaps to only remember these faces that would be cut
			//these are the faces that should be trimmed to avoid collision

			return 0;
		}

		//find optimal cuts of the net to avoid collision
		//return cost of cutting the net
		float BloomingUnfolding::Net::optimalcuts(Unfolder *unfolder, const vector<set<uint>>& overlaps)
		{
			model * m = unfolder->getModel();
			BitVector killed(m->t_size); //record what face is removed

			vector<float> face_areas(m->t_size, 0); //
			compute_area_in_subnet(m, this->root, face_areas);
			list<int> Q; //a queue
			Q.push_back(unfolder->getBaseFaceID());
			float cost = optimalcuts(unfolder, Q, this->code, killed, face_areas, overlaps);

			//this is a temporary patch....for the diagonal edges that are cut due to intersection....
			//diagonal edges are not supposed to be cut so we need to find a non-diagonal edge (likely its parent)
			//Note: this should be handled in more general ways
			for (int i = 0; i < m->e_size; i++) {
				if (this->code[i] != 0) continue;
				edge& e = m->edges[i];
				if (e.type == 'd') {
					//find the parent....
					uint f1 = e.fid.front();
					uint f2 = e.fid.back();
					uint f1p = parent[f1];
					uint p = (f2 == f1p) ? f2 : f1;
					uint pe= m->getEdgeIdByFids(p,parent[p]);
					if (pe == -1) {
						cerr << "! Error: BloomingUnfolding::Net::optimalcuts" << endl;
						exit(1);
					}
					cout << "cut " << pe << " because " << i << " is cut and is a diagonal edge" << endl;
					code.off(pe);
				}
			}

			rebuild(unfolder, this->code);//rebuild this net from code

			return cost;
		}

		//find optimal cuts of the net
		//return cost of cutting the net
		//the resulting net is stored in B
		float BloomingUnfolding::Net::optimalcuts
		(Unfolder *unfolder, list<int>& Q, BitVector& B, BitVector& killed,
		 const vector<float>& face_areas, const vector<set<uint>>& _O)
		{
			if (Q.empty()) return 0; ///nothing to see here

			model * m = unfolder->getModel();
			int f = Q.front();
			Q.pop_front();

			if (killed[f]) {
				float cost= optimalcuts(unfolder, Q, B, killed, face_areas, _O);
				return cost;
			}

			//original overlaps
			const set<uint>& overlaps_orig = _O[f];

			list<uint> overlaps; //updated overlaps
			for (uint o : overlaps_orig) {
				if (!killed[o]) overlaps.push_back(o);
			}

			if (overlaps.empty()) { //collision free face
				for (uint kid : this->kids[f]) {
					int eid = m->getEdgeIdByFids(f, kid);
					const edge& e = m->edges[eid];
					if (e.type == 'd' || B[eid]) Q.push_back(kid);
				}
				float cost= optimalcuts(unfolder, Q, B, killed, face_areas, _O);
				return cost;
			}

			//case 1, remove the colliding face
			float cost1 = 0;
			BitVector B1 = B;
			BitVector K1 = killed;
			{
				list<int> Q1 = Q;

//				cout<<"at F="<<f<<" we kill f=";
				for (uint id : overlaps) { //remove faces
					uint eid = m->getEdgeIdByFids(id, parent[id]);
					B1.off(eid);
					killDecendents(id,K1);  //mark the face and its desendents as killed so we don't hanlde them

					auto ftype = (BloomingUnfolding::FACE_TYPE)m->tris[id].cluster_id;
					cost1 += g_face_cost.at(ftype)*face_areas[id];
				}
				//cout<<endl;

				for (uint kid : this->kids[f]) {
					int eid = m->getEdgeIdByFids(f, kid);
					const edge& e = m->edges[eid];
					if (e.type == 'd' || B[eid]) Q1.push_back(kid);
				}
				if (!Q1.empty())
					cost1 += optimalcuts(unfolder, Q1, B1, K1, face_areas, _O);
			}

			//case 2, remove f
			float cost2 = 0;
			BitVector B2 = B;
			BitVector K2 = killed;
			{
				uint eid = m->getEdgeIdByFids(f, parent[f]);
				B2.off(eid);
				killDecendents(f, K2);

				cost2 = optimalcuts(unfolder, Q, B2, K2, face_areas, _O);
				auto ftype = (BloomingUnfolding::FACE_TYPE)m->tris[f].cluster_id;
				cost2 += g_face_cost.at(ftype)*face_areas[f];
			}

			//compare
			if (cost1 < cost2) {
				B = B1;
				return cost1;
			}
			//else
			B = B2;
			return cost2;
		}

		//compute a list of areas that is in the subnet of a give face
		//this value include the area of the face itself
		//pre consition: vector<float>& areas must contain m->t_size zeros when this function is first
		//called
		void BloomingUnfolding::Net::compute_area_in_subnet(model *m, int root, vector<float>& areas)
		{
			if (areas.size() != m->t_size || areas[root] != 0)
				areas = vector<float>(m->t_size, 0);

			triangle& tri = m->tris[root];
			//areas[root] = tri.area;

			auto type=(BloomingUnfolding::FACE_TYPE)tri.cluster_id;
			auto cost=BloomingUnfolding::g_face_cost.at(type);
			areas[root] = cost*tri.area;

			for (int kid : this->kids[root])
			{
				compute_area_in_subnet(m, kid, areas);
				areas[root] += areas[kid];
			}//edn for i
		}

		//compute a list of cost for flipping an edge
		//calls compute_edge_cost(model *m, int face, vector<float>& costs)
		void BloomingUnfolding::Net::compute_edge_cost(model *m, vector<float>& costs)
		{
			costs=vector<float>(m->e_size,0);
			for (int kid : this->kids[this->root])
			{
				compute_edge_cost(m, kid, costs);
			}//edn for i
		}

		//compute a list of cost for flipping an edge
		//the cost is based on the sum of its faces in the subtree
		//the edge is defined as the edge between "face" and its parent
		float BloomingUnfolding::Net::compute_edge_cost(model *m, int face, vector<float>& costs)
		{
			assert(face!=this->root);
			int eid=m->getEdgeIdByFids(face, this->parent[face]);
			const triangle& t=m->tris[face];
			auto type=(BloomingUnfolding::FACE_TYPE)t.cluster_id;
			auto cost=BloomingUnfolding::g_face_cost.at(type);

			costs[eid] = cost*t.area;
			for (int kid : this->kids[face])
			{
				costs[eid] += compute_edge_cost(m, kid, costs);
			}//edn for i

			//cout<<"costs["<<eid<<"]="<<costs[eid]<<endl;
			return costs[eid];
		}



		//compute the type of each face in the model that is connected to the root
		//in the given net
		void BloomingUnfolding::Net::compute_face_types(Unfolder * unfolder) const
		{
			model * m = unfolder->getModel();

			//collect all keep faces connected to the root
			//they will be marked as BASE_FACE
			int fbase = unfolder->getBaseFaceID();
			m->tris[fbase].cluster_id = BASE_FACE;
			list<int> open;
			open.push_back(fbase);
			while(!open.empty()){
				int fid=open.front();
				open.pop_front();
				m->tris[fid].cluster_id = BASE_FACE;
				for(int kid : kids[fid]){
					auto ftype = (BloomingUnfolding::FACE_TYPE)m->tris[kid].cluster_id;
					if(ftype==KEEP_FACE) open.push_back(kid);
				}
			}//end whiel

			//mark every BASE_FACE a KEEP_FACE
			for(int i=0;i<m->t_size;i++){
				auto ftype = (BloomingUnfolding::FACE_TYPE)m->tris[i].cluster_id;
				if(ftype==BASE_FACE){
					m->tris[i].cluster_id = KEEP_FACE;
				}
				else m->tris[i].cluster_id=TOSS_FACE;
			}
			//mark the true base face
			m->tris[fbase].cluster_id = BASE_FACE;
		}

		void BloomingUnfolding::Net::killDecendents(uint fid, BitVector& killed) {
			killed.on(fid);
			for (uint kid : this->kids[fid])  killDecendents(kid,killed);
		}

		//
		//
		// FlowerBloomingSplitter (Not USED)
		//
		//


		FlowerBloomingSplitter::FlowerBloomingSplitter(model* m, int base) {
			m_base = base; m_model = m;
			assert(m_model);
		}


		Vector3d FlowerBloomingSplitter::genRandomUnitVector(const Config& config)
		{
			//cout<<"-m_model->tris[m_base].n="<<-m_model->tris[m_base].n<<endl;
			return -m_model->tris[m_base].n;
		}

		//get a list of  faces coplanar to the given face
		//the return list include the face itself
		//Note: this should actually be model to model.h
		inline void getCoplanarFaces(model * m, int face, list<int>& coplanars)
		{
			list<int> open;
			set<int> closed;
			open.push_back(face);

			//get the faces coplanar to the base
			while (!open.empty()) {
				int fid = open.front();
				open.pop_front();
				triangle & f = m->tris[fid];
				closed.insert(fid);

				for (int i = 0; i < 3; i++) {
					edge& e = m->edges[f.e[i]];
					//const BitVector& code=net.encode();
					if (e.type != 'd') continue;
					int of = e.otherf(fid);
					if (closed.find(of) != closed.end()) continue; //this can be slow, but we should not have too many faces in closed
					if (of != -1) open.push_back(of);
				}//end for i
			}//end while

			coplanars.insert(coplanars.end(), closed.begin(), closed.end());
		}

	}
}
