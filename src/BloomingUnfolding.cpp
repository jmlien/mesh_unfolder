/*
* ClusterUnfolding.cpp
*
*  Created on: Dec 01, 2015
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
		const map<BloomingUnfolding::FACE_TYPE, float> BloomingUnfolding::g_face_cost =
		{ {BloomingUnfolding::BASE_FACE, 10000.0f},
		 {BloomingUnfolding::KEEP_FACE, 1.0f},
					 {BloomingUnfolding::TOSS_FACE, 0.001f} };

		//
		list<int> neighbor_edges(const BitVector& code, model * m, int eid);

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
					if (e.type != 'd') continue;
					int of = e.otherf(fid);
					if (closed.find(of) != closed.end()) continue; //this can be slow, but we should not have too many faces in closed
					if (of != -1) open.push_back(of);
				}//end for i
			}//end while

			coplanars.insert(coplanars.end(), closed.begin(), closed.end());
		}
/*
		const int BitVector::_S = sizeof(int) * 8; //number of bits for int

		BitVector::BitVector(const vector<bool>& values)
			:_N((int)values.size())
		{
			data = vector<unsigned int>((int)ceil(_N*1.0f / _S), 0);

			for (int i = 0; i < _N; i++) {
				//cout<<"trying i="<<i<<endl;
				if (values[i]) on(i);
				// if((*this)[129]){
				//   cout<<"129 is on why? "<<i<<endl;
				//   exit(1);
				// }
			}
		}

		BitVector::BitVector(int size)
			:_N(size)
		{
			data = vector<unsigned int>((int)ceil(_N*1.0f / _S), 0);
		}

		void BitVector::on(int i) {
			assert(i < _N);
			data[i / _S] |= (1 << (i%_S));
		}

		void BitVector::off(int i) {
			assert(i < _N);
			data[i / _S] &= ~(1 << (i%_S));
		}

		bool BitVector::operator[](int i) const {
			assert(i < _N);
			return data[i / _S] & (1 << (i%_S));
		}

		bool BitVector::operator<(const BitVector& other) const {
			int size = data.size();
			for (int i = size - 1; i >= 0; i--) {
				if (data[i] != other.data[i]) return data[i] < other.data[i];
			}
			return false; //they are the same
		}

		ostream& operator<<(ostream& out, const BitVector& v) {
			out << "(";
			for (int i = 0; i < v._N; i++) out << v[i] << ((i < v._N - 1) ? "," : ")");
			return out;
		}

		bool BitVector::operator==(const BitVector& other) const {
			if (other._N != _N) return false;
			int size = data.size();
			for (int i = 0; i < size; i++) if (data[i] != other.data[i]) return false;
			return true;
		}

		void BitVector::tovector(vector<bool>& values) const {
			values.resize(_N); //ensure that there is enough space
			for (int i = 0; i < _N; i++) values[i] = (*this)[i];
		}
*/
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

			m_blooming_method = FLOWER;
			m_splitter = NULL;
			m_blooming_base_fase = mycfg.baseface;
			m_runs = mycfg.max_retries;
			m_start_time = 0;
			this->m_unfolder = unfolder;
			this->m_preview_only = false;
			this->m_no_trim=false;
		}

		BloomingUnfolding::~BloomingUnfolding()
		{

		}

		bool BloomingUnfolding::setup(const string& method)
		{

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

			cout<<"(";
			if (method.find("preview") != string::npos)
			{
				this->m_preview_only = true; //enable preview only mode
				cout << "\"preview\" the initial unfolding only;" << endl;
			}
			else if (method.find("nocut") != string::npos)
			{
				this->m_no_trim = true; //avoid trimming the model
				cout << "\"preview\" A* unfolding only;" << endl;
			}
			else{
				cout << "add \"preview\" to see the initial unfolding;";
				cout << "add \"nocut\" to see the overlapping faces;";
			}
			cout<<")"<<endl;

			assert(m_splitter);
			if (m_blooming_base_fase < 0 || m_blooming_base_fase >= m_unfolder->getModel()->t_size) {
				cerr << "! Error: BloomingUnfolding::setup: base face id=" << m_blooming_base_fase << " is out of bound" << endl;
				return false;
			}
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

		void BloomingUnfolding::run()
		{
			const auto config = this->m_unfolder->getConfig();

			m_start_time = clock();

			bool successful = false;

			auto mesh = this->m_unfolder->getModel();

			m_splitter->measure(mesh);

			Config mycfg = m_unfolder->getConfig();
			vector<float> best_weights;
			int min_overlap_count = INT_MAX;
			for (int r = 0; r < this->m_runs; r++)
			{
				//cfg.user_vector = rand_dirs[r];
				vector<float> weights = m_splitter->assignWeights(mesh, mycfg);
				auto count = m_unfolder->buildFromWeights(weights);
				if (count < min_overlap_count) {
					best_weights = weights;
					min_overlap_count = count;
					cerr << "- iter = " << r << ", total overlaps = " << count << "\r" << flush;
					if (count == 0) break; //done, found a net
				}
			}
			cerr<<"\n"<<flush;

			//finalize
			m_unfolder->buildFromWeights(best_weights);
			if (this->m_preview_only) return; //preview the unfolding without optimization

#if 1
			//
			//get creases from the current unfolding
			//we assume that this is the best blooming unfolding
			//
			const set<uint>& creases = m_unfolder->getFoldEdges();
			vector<bool> ic(mesh->e_size, 0);
			for (uint c : creases) ic[c] = true;

			//create a net from the creases
			Net start(mesh, m_unfolder->getBaseFaceID(), BitVector(ic));
			AStar_Helper_Mixed mixed(start, m_unfolder, 1);

			Net optmized_net = AStar(mixed);
#else
			Net optmized_net = AStar();
#endif

			m_unfolder->buildFromWeights(optmized_net.encode());
/*
			auto& overlaps = m_unfolder->getOverlppingFacePairs();
			for(int i=0;i<mesh->t_size;i++){
				cout<<"face["<<i<<"] intersects ";
				for(auto& o : overlaps[i]){
					cout<<o<<";";
				}
				cout<<endl;
			}
*/

			//remove faces that are overlapping....
			if(!this->m_no_trim){ //this is default
				optmized_net.optimalcuts(m_unfolder);
				m_unfolder = build_unfolder_from_net(m_unfolder, optmized_net);
				model * subm=m_unfolder->getModel();
				stringstream ss;
				ss << "./" << mesh->name << "_bloom_" << subm->t_size << ".obj";
				subm->saveObj(ss.str());
				cout<<"- Blooming unfolding saved the trimmed model to "<<ss.str()<<endl;
			}

		}// END BloomingUnfolding::run()


		BloomingUnfolding::Net BloomingUnfolding::AStar(AStar_Helper& help)
		{
			//build the first net from unfolder
			model * m = m_unfolder->getModel();
			BloomingUnfolding::g_dist_level = 1; //start with level 1: full distance
			const int fbase = m_unfolder->getBaseFaceID(); //base face of a net, same of all nets

			//compute addition data to help guiding A* search
			compute_face_types(m_unfolder, m->tris[fbase].n);
			vector<float> ecosts; //cost of flipping a root edge
			help.root.compute_edge_cost(m, ecosts);

			//start A* seach from the start net
			vector<Net> open;
			set<BitVector> visted;

			Net best_net = help.root;
			visted.insert(help.root.encode());
			open.push_back(help.root);

			const int max_failed_attempts = 1000;
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
				else if (help.heuristics(net) < help.heuristics(best_net))
				{
						best_net = net;
						cout << "- Best net overlaps=" << help.heuristics(net)
							<< " dist=" << help.dist(net,ecosts) << " open size=" << open.size() << endl;
						failed_attempts = 0; //reset
				}
				else {//failed to find
					failed_attempts++;
					if (failed_attempts % (max_failed_attempts/10) == 0)
						cout << "! Failed attempt count: " << failed_attempts << "/" << max_failed_attempts << endl;

					if (failed_attempts > max_failed_attempts) {

						if (BloomingUnfolding::g_dist_level == 3) { //still failed as max level?
							cerr << "! Max level reached. Report the best unfolding." << endl;
							return best_net;
						}

						cout << "- Max failed attempts reached (" << max_failed_attempts << "); change distance matrics to level " << BloomingUnfolding::g_dist_level + 1
							<< " open size=" << open.size() << endl;

						//OK, move to the next level of distance metric
						BloomingUnfolding::g_dist_level++;
						BloomingUnfolding::g_recompute_dist = true;
						for (Net& net : open) help.dist(net,ecosts); //update the distance
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
			Net start(m, fbase, BitVector(ic));

			//compute addition data to help guiding A* search
			vector<float> ecosts; //cost of flipping a root edge
			compute_face_types(m_unfolder, m->tris[fbase].n);
			start.compute_edge_cost(m, ecosts);

			//start A* seach from the start net
			vector<Net> open;
			Net best_net = start;

			//map<BitVector, float> gscore;
			//map<BitVector, float> fscore;
			set<BitVector> visted;
			//set<BitVector> closed;
			//unordered_set<BitVector,hash_BitVector> visted;
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

						cout << "- Max failed attempts reached (" << max_failed_attempts << "); change distance matrics to level " << BloomingUnfolding::g_dist_level + 1
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
				list<Net> neighbors = net.neighbors(m, fbase);
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

		//build a model from the bit vector
		//this likely is used when some faces of the input net is trimmed by the ned
		//the bitvector may be updated to reflect the new model
		model * BloomingUnfolding::build_model_from_net(model * m, int base, BitVector& B)
		{
			masc::obj::objModel data;

			//we will need a list of faces first, if the number if the same as m->t_size
			//then we are done.
			list<uint> fids;
			{//this should be moved to a function actually
				list<int> open;
				open.push_back(base);
				vector<bool> visited(m->t_size, false);


				while (!open.empty()) {
					int n = open.front();
					open.pop_front();
					fids.push_back(n); //remember all the faces we visited
					visited[n] = true;

					for (int i = 0; i < 3; i++) {
						int eid = m->tris[n].e[i];

						edge& e = m->edges[eid];
						if (!B[eid] && e.type != 'd') continue; //not connected

						int of = e.otherf(n);
						if (visited[of]) continue; //visited
						open.push_back(of);
					}//edn for
				}//end while
			}


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
			if (subm->build(data, true) == false)
			{
				cerr << "! Error: Failed to build a model from net" << endl;
				return NULL;
			}

			//register the new fid to the old fid
			auto fid_it = fids.begin();
			for (int fid = 0; fid < subm->t_size; fid++)
			{
				subm->tris[fid].source_fid = *fid_it;
				++fid_it;
			}


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

		//compute the type of each face in the model using the given viewing dirsection
		//store face type in data "cluster_id"
		void BloomingUnfolding::compute_face_types
		(Unfolder * unfolder, const Vector3d& dir)
		{
			model * m = unfolder->getModel();

			for (int i = 0; i < m->t_size; i++) {
				triangle& tri = m->tris[i];
				if (i == unfolder->getBaseFaceID()) { tri.cluster_id = BASE_FACE; }
				else {
					//we abuse cluster_id here
					tri.cluster_id = (dir*tri.n >= 0) ? KEEP_FACE : TOSS_FACE;
				}//end if
			}//end for i
		}

//
//
//
//
//

//-------------------
list<BloomingUnfolding::Net>
BloomingUnfolding::AStar_Helper_Mixed::neighbors(BloomingUnfolding::Net& net)
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
			neis.push_back(Net(m, fbase, ncode));
		}
	}//end for i

	return neis; //netghbors of this net
}

float BloomingUnfolding::AStar_Helper_Mixed::dist(BloomingUnfolding::Net& net, const vector<float>& ecosts)
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

float BloomingUnfolding::AStar_Helper_Mixed::heuristics(BloomingUnfolding::Net& net)
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
list<BloomingUnfolding::Net> BloomingUnfolding::AStar_Helper_Keep::neighbors(BloomingUnfolding::Net& net)
{
	list<Net> neis;
	return neis; //netghbors of this net
}

float BloomingUnfolding::AStar_Helper_Keep::dist(BloomingUnfolding::Net& net, const vector<float>& ecosts)
{

	return 0;
}

float BloomingUnfolding::AStar_Helper_Keep::heuristics(BloomingUnfolding::Net& net)
{

	return 0;
}

//-------------------
list<BloomingUnfolding::Net> BloomingUnfolding::AStar_Helper_Toss::neighbors(BloomingUnfolding::Net& net)
{
	list<Net> neis;
	return neis; //netghbors of this net
}

float BloomingUnfolding::AStar_Helper_Toss::dist(BloomingUnfolding::Net& net, const vector<float>& ecosts)
{
	return 0;
}

float BloomingUnfolding::AStar_Helper_Toss::heuristics(BloomingUnfolding::Net& net)
{
	return 0;
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
			this->depth = other.depth;
			this->leaves = other.leaves;
			this->parent = other.parent;
			this->kids = other.kids;
		}

		// BloomingUnfolding::Net::Net(model * m, int root, const vector<bool>& value) //number of edges in this net
		// 	:BloomingUnfolding::Net::Net(m,root,BitVector(value))
		// {}

		BloomingUnfolding::Net::Net(model * m, int root, const BitVector& value) //number of edges in this net
			:code(value)
		{
			//init data
			this->root = root;
			g = h = FLT_MAX;
			this->n = code.size();
			analyze(m, code);
		}

		//get the neighbors of this net
		list<BloomingUnfolding::Net> BloomingUnfolding::Net::neighbors(model * m, int fbase)
		{
			list<Net> neis;
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

					neis.push_back(Net(m, fbase, ncode));
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
		void BloomingUnfolding::Net::analyze(model *m, const BitVector& value)
		{
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
		void BloomingUnfolding::Net::rebuild(model * m, const BitVector& B)
		{
			//init data
			this->code = B;
			this->depth = 0;
			this->leaves = 0;
			this->parent.clear();
			this->kids.clear();

			g = h = FLT_MAX;
			this->n = (int)this->code.size();
			analyze(m, this->code);
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
			/*for (uint i = 0; i < m->t_size; i++) {//for each face
				if (overlaps[i].empty()) continue; //nothing here, move on
				const triangle& ti = m->tris[i];
				auto type_i=(BloomingUnfolding::FACE_TYPE)ti.cluster_id;
				auto cost_i=BloomingUnfolding::g_face_cost.at(type_i);
				cost += cost_i;
			}
			*/

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

		//find optimal cuts of the net
		//return cost of cutting the net
		//the resulting net is stored in B
		float BloomingUnfolding::Net::optimalcuts(Unfolder *unfolder)
		{
			model * m = unfolder->getModel();
			BitVector killed(m->t_size); //record what face is removed

			vector<float> face_areas(m->t_size, 0); //
			compute_area_in_subnet(m, this->root, face_areas);
			list<int> Q; //a queue
			Q.push_back(unfolder->getBaseFaceID());
			float cost = optimalcuts(unfolder, Q, this->code, killed, face_areas);

			rebuild(m, this->code);//rebuild this net from code
			return cost;
		}

		//find optimal cuts of the net
		//return cost of cutting the net
		//the resulting net is stored in B
		float BloomingUnfolding::Net::optimalcuts(Unfolder *unfolder, list<int>& Q, BitVector& B, BitVector& killed, const vector<float>& face_areas)
		{
			if (Q.empty()) return 0; ///nothing to see here

			model * m = unfolder->getModel();
			int f = Q.front();
			Q.pop_front();



			const set<uint>& overlaps = unfolder->getOverlppingFacePairs()[f];


			if (killed[f]) {
				return optimalcuts(unfolder, Q, B, killed, face_areas);
			}

			if (overlaps.empty()) { //collision free face
				for (uint kid : this->kids[f]) {
					int eid = m->getEdgeIdByFids(f, kid);
					const edge& e = m->edges[eid];
					if (e.type == 'd' || B[eid]) Q.push_back(kid);
				}
				return optimalcuts(unfolder, Q, B, killed, face_areas);
			}

			//case 1, remove the colliding face
			float cost1 = 0;
			BitVector B1 = B;
			{
				BitVector K1 = killed;
				list<int> Q1 = Q;

//				cout<<"at F="<<f<<" we kill f=";
				for (uint id : overlaps) { //remove faces
					uint eid = m->getEdgeIdByFids(id, parent[id]);
					B1.off(eid);
					K1.on(id); //mark the face as killed so we don't hanlde it
					//cout<<id<<";";
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
					cost1 += optimalcuts(unfolder, Q1, B1, K1, face_areas);
			}

			//case 2, remove f
			float cost2 = 0;
			BitVector B2 = B;
			{
				uint eid = m->getEdgeIdByFids(f, parent[f]);
				B2.off(eid);
				cost2 = optimalcuts(unfolder, Q, B2, killed, face_areas);
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

		//
		//
		// FlowerBloomingSplitter
		//
		//


		FlowerBloomingSplitter::FlowerBloomingSplitter(model* m, int base) {
			m_base = base; m_model = m;
			assert(m_model);
		}


		Vector3d FlowerBloomingSplitter::genRandomUnitVector(const Config& config)
		{
			return -m_model->tris[m_base].n;
		}

		StarBloomingSplitter::StarBloomingSplitter(model* m, int base)
			:FlowerBloomingSplitter(m, base)
		{
			//find vertices at the base
			list<int> coplanars;
			getCoplanarFaces(m, base, coplanars);
			for (int fid : coplanars) {
				for (int i = 0; i < 3; i++) m_base_vids.insert(m->tris[fid].v[i]);
			}

			//now find the cuts by building the dijstra tree
			vector< pair<float, int> > open;
			for (int vid : m_base_vids) {
				vertex& v = m_model->vertices[vid];
				v.source_vid = vid;
				v.cut_src_id = vid;
				v.score = 0;
				open.push_back(make_pair(0, vid));
			}
			make_heap(open.begin(), open.end());
			//run dijstra from vertices of the base face
			while (!open.empty()) {

				pair<float, int> o = open.front();
				pop_heap(open.begin(), open.end());
				open.pop_back();
				int vid = o.second;
				vertex& v = m_model->vertices[vid];
				if (o.first > v.score) continue; //out of date data

				//spread to the neighbors
				for (int eid : v.m_e) {
					edge& e = m_model->edges[eid];
					int nid = e.otherv(vid);
					vertex& n = m_model->vertices[nid];
					bool updated = false;

					if (n.cut_src_id == UINT_MAX) updated = true;//not visited
					else if (n.score > e.length + v.score) updated = true;//visited, but found a shorter path

					if (updated) {
						n.score = e.length + v.score;
						n.cut_src_id = vid; //parent
						n.source_vid = v.source_vid; //root
						open.push_back(make_pair(n.score, nid));
						push_heap(open.begin(), open.end());
					}

				}//end for nei
			}//end while

			//go through all vertices and collect local maxium
			list<int> local_max_vids;

			for (int i = 0; i < m_model->v_size; i++) {
				vertex& v = m_model->vertices[i];
				cout << "v[" << i << "].score=" << v.score << endl;
				//check if v is local max
				bool local_max = true;
				for (int eid : v.m_e) {
					vertex& n = m_model->vertices[m_model->edges[eid].otherv(i)];
					if (n.score > v.score) {
						local_max = false;
						break;
					}
				}//end for eid

				if (local_max)
					local_max_vids.push_back(i);
			}//end for i

			//trace back to the base face
			m_cut_edges.clear();
			for (int vid : local_max_vids) {
				cout << "local max=" << vid << endl;
				trace(vid, m_cut_edges);
				vertex& v = m_model->vertices[vid];
				//check if its neighor can lead to other source in the same
				//distance
				for (int eid : v.m_e) {
					int nid = m_model->edges[eid].otherv(vid);
					if (nid == v.cut_src_id) continue; //parent
					vertex& n = m_model->vertices[nid];
					if (n.score + m_model->edges[eid].length == v.score && v.source_vid != n.source_vid) {
						v.cut_src_id = nid;
						trace(vid, m_cut_edges);//trace again with this new parent
					}
				}//end for eid
			}//end of vid

			for (int eid : m_cut_edges) {
				edge& e = m_model->edges[eid];
				cout << "cut (" << e.fid.front() << "," << e.fid.back() << ")" << endl;
			}

			//reset
			for (int i = 0; i < m_model->v_size; i++) {
				vertex& v = m_model->vertices[i];
				v.score = 0;
				v.cut_src_id = UINT_MAX;
				v.source_vid = UINT_MAX;
			}//end for i

		}

		//trace a list of edges
		void StarBloomingSplitter::trace(int vid, unordered_set<int>& eids)
		{
			while (true) {
				vertex& v = m_model->vertices[vid];
				int p = v.cut_src_id;
				if (p == vid) break; //root here
				eids.insert(m_model->getEdgeId(vid, p));
				vid = p;
			}
		}

		void StarBloomingSplitter::assignWeightsImpl(model* m, vector<float>& weights, const Config& config)
		{
			FlowerBloomingSplitter::assignWeightsImpl(m, weights, config);
			for (int eid : m_cut_edges) weights[eid] = 10;
		}

	}
}
