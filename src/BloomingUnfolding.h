/*
* BloomingUnfolding.h
*
*  Created on: April 22, 2020
*      Author: jmlien
*/

#pragma once

#include <string>
#include <vector>

#include "unfolder.h"
#include "Splitter.h"
#include "util/BitVector.h"

namespace masc {
	namespace unfolding {

		using namespace std;
/*
		//a long stream of bits
		class BitVector {
		public:

			BitVector(const vector<bool>& values);
			BitVector(int size);//initialize all to 0s

			int size() const { return _N; }
			void on(int i); //turn on i-th bit
			void off(int i); //turn off i-th bit

			//operator overloading
			bool operator[](int i) const;
			bool operator<(const BitVector& other) const; //comparator
			friend ostream& operator<<(ostream& out, const BitVector& v);
			bool operator==(const BitVector& other) const;

			//convert BitVector to a vector of bool
			void tovector(vector<bool>& values) const;

		protected:
			vector<unsigned int> data;
			const static int _S; //this is sizeof(unsigned int)
			int _N; //number of bits stored
			friend struct hash_BitVector;
		};

		struct hash_BitVector {
			//for hashing
			size_t operator()(const BitVector &v) const
			{
				size_t code = 0x2D2816FE;
				for (int i : v.data) return code * 31 + i; //hash<int>()(i);
				return code;
			}
		};
*/
		class BloomingUnfolding
		{
		public:

			BloomingUnfolding(Unfolder* unfolder);
			virtual ~BloomingUnfolding();

			bool setup(const string& method);

			//////////////////////////////////////////
			// run n unfoldings and cluster
			//////////////////////////////////////////

			virtual void run();

			//////////////////////////////////////////

			virtual void print(ostream& out) const;

			//////////////////////////////////////////
			// access
			//////////////////////////////////////////
			Unfolder* getUnfolder() { return m_unfolder; }

		protected:

			// will be called once after setup,
			// can be used to override parameters
			void virtual init();


			//perform blooming unfolding
			bool blooming_unfold(Unfolder * unfolder);

			class Net //represented by a binary variable
			{
			public:

				Net(const Net& other); //copy constructor
				Net(model * m, int root, const BitVector& value); //build from value

				//get the neighbors of this net
				//requires the model and the base face id
				list<Net> neighbors(model * m, int fbase);

				//the distance from the root
				float dist(model * m, const Net& root, vector<float>& ecosts);

				//distance to go heuristics
				float heuristics(Unfolder * unfolder);

				//convert this net to a hash code
				const BitVector& encode() const;

				//comparator
				bool operator<(const Net& other) const;

				//given a net find an optimal way to remove overlapping faces
				float optimalcuts(Unfolder *unfolder);

				//compute a list of cost for flipping an edge
				//calls compute_edge_cost(model *m, int face, vector<float>& costs)
				void compute_edge_cost(model *m, vector<float>& costs);


				//compute the type of each face in the model that is connected to the root
				//in the given net
				void compute_face_types(Unfolder * unfolder) const;

				//access
				float & G(){ return g; }
				float & H(){ return h; }
				float getDepth() const { return depth; }
				float getLeafSize() const { return leaves; }
				//const list<int>& Kids(int fid) const { return kids[fid]; }
				//int Parent(int fid) const { return parent[fid]; }

			private:

				float g, h; //g, h scores used in A*, f=g+h
				float depth, leaves; //number of leaves and depth of the net
				int n; //=code.size(), this should be the same for all nets
				BitVector code;

				//tree
				int root; //root id, this should be the same for all nets
				vector<int> parent;      //parent[i] is the parent of face i
				vector<list<int> > kids; //kids[i] is the kids of face i

			protected:

				//a naive method to
				//compute a list of areas that is in the subnet of a give face
				//this value include the area of the face itself
				//pre consition: vector<float>& areas moptimalcutsust contain m->t_size zeros when this function is first
				//called
				//float overlapping_cost(Unfolder *unfolder);

				//find optimal cuts of the net
				//return cost of cutting the net
				//the resulting net is stored in B
				float optimalcuts(Unfolder *unfolder, list<int>& Q, BitVector& B,
					                BitVector& killed,
													const vector<float>& areas);

				//compute a list of areas that is in the subnet of a give face
				//this value include the area of the face itself
				void compute_area_in_subnet(model *m, int root, vector<float>& areas);

				//compute a list of cost for flipping an edge
				//the cost is based on the sum of its faces in the subtree
				//the edge is defined as the edge between "face" and its parent
				float compute_edge_cost(model *m, int face, vector<float>& costs);

				//get a list of edges if the eid is connected in n
				//list<int> neighbor_edges(model * m, int eid);

				//compute some statistics of this net, directly from code
				//count # of leaves  & distance from the base to the farthest leaf
				void analyze(model *m, const BitVector& value);

				//rebuild this net from the give code
				void rebuild(model * m, const BitVector& B);
			};

			struct AStar_Helper
			{
				AStar_Helper(const Net& net, Unfolder * uf, float t=0, int ml=3)
				: root(net){ termination=t; unfolder=uf; max_level=ml;}
				virtual list<Net> neighbors(Net& net)=0;
				virtual float dist(Net& net, const vector<float>& ecosts)=0;
				virtual float heuristics(Net& net)=0;
				Unfolder * unfolder;
				Net root;
				float termination;
				int max_level;
			};

			//optimize the net with A*
			Net AStar();//this is the first AStar, which is the same as AStar mix

			Net AStarMix();//same as AStar

			Net AStar2Steps(); //try to unfolde the keep faces and then toss faces

			//optimize the net with A*
			Net AStar(AStar_Helper & help);

			//build an unfolder that has only the KEEP_FACEs
			Unfolder * build_keep_unfolder(Unfolder *unfolder, Net& net);

			//build a new unfolder from the net, this is used most likely that
			//some faces should be removed to ensure that there is no overlapping
			Unfolder * build_unfolder_from_net(Unfolder *unfolder, Net& net);

			//build a model from the bit vector
			//this likely is used when some faces of the input net is trimmed by the ned
			//the bitvector may be updated to reflect the new model
			model * build_model_from_net(model * m, int base, BitVector& B);

			//compute the type of each face in the model using the given viewing dirsection
			//the computed information is stored in model.triangle.cluster_id
			void compute_face_types(Unfolder * unfolder, const Vector3d& dir);

			//create a subset of m using the fids
			model * build_submodel(model * m, list<uint> fids);

			//get a list of faces that are connected to the base face according to the BitVector
			//faces of the given toss_face_type are ignored
			list<int> getCC(model * m, int base, const BitVector& B, int toss_face_type=-1);

			struct AStar_Helper_Mixed : public AStar_Helper
			{

				AStar_Helper_Mixed(const Net& net, Unfolder * uf, float t=0, int ml=3):AStar_Helper(net,uf,t,ml){}
				list<Net> neighbors(Net& net) override;
				float dist(Net& net, const vector<float>& ecosts) override;
				float heuristics(Net& net) override;
			};

			struct AStar_Helper_Keep : public AStar_Helper
			{
				AStar_Helper_Keep(const Net& net, Unfolder * uf, float t=0, int ml=3):AStar_Helper(net,uf,t,ml){}
				list<Net> neighbors(Net& net) override;
				float dist(Net& net, const vector<float>& ecosts) override;
				float heuristics(Net& net) override;
				bool isKeep(int eid); //is the edge incident to keep faces
			};

			struct AStar_Helper_Toss : public AStar_Helper
			{
				AStar_Helper_Toss(const Net& net, Unfolder * uf, float t=0, int ml=3):AStar_Helper(net,uf,t,ml){}
				list<Net> neighbors(Net& net) override;
				float dist(Net& net, const vector<float>& ecosts) override;
				float heuristics(Net& net) override;
				bool isKeep(int eid); //is the edge incident to keep faces
				//bool Keep_KID(BloomingUnfolding::Net& net, int f);
				//bool Keep_SUBNET(BloomingUnfolding::Net& net, int e);
			};

		private:

			//
			//enum BLOOMING_METHOD { FLOWER, STAR };
			//BLOOMING_METHOD m_blooming_method;

			//Splitter used in this unfolding
			//Splitter * m_splitter;

			//blooming base face
			//int m_blooming_base_fase; //user provided

			uint m_runs; //number of trys

			long m_start_time;

			Unfolder* m_unfolder; //a reference to the unfolder provided

			bool m_preview_only; //no A* optimization

			bool m_no_trim; //do not trim after A*

			static bool g_recompute_dist;
			static short g_dist_level; //1~3, level (1) dist with all matrics, (2) dist with bit flips (3) no dist

			//types of faces
			//this information is stored in model.triangle.cluster_id
			enum FACE_TYPE {BASE_FACE, KEEP_FACE, TOSS_FACE};

			//cost for removing each type of face
			static map<FACE_TYPE, float> g_face_cost; //cost for removing a face from the net

			//
			static float overlapping_cost(Unfolder *unfolder);
		};


		class FlowerBloomingSplitter : public SteepestEdgeSplitter
		{
		public:

			FlowerBloomingSplitter(model* m, int base);
			virtual ~FlowerBloomingSplitter() {};
			string getName() override { return "FlowerBloomingSplitter"; }
			Vector3d genRandomUnitVector(const Config& config) override;

		protected:

			//same as SteepestEdgeSplitter
			//void assignWeightsImpl(model* m, vector<float>& weights, const Config& config) override;
			int m_base;
			model * m_model;
		};

/*
		class StarBloomingSplitter : public FlowerBloomingSplitter
		{
		public:
			StarBloomingSplitter(model* m, int base);
			virtual ~StarBloomingSplitter() {};
			string getName() override { return "StarBloomingSplitter"; }
		protected:

			set<int> m_base_vids; //vids from the base faces

			//trace a path that links to the root vertex
			//store the trace in eids
			void trace(int vid, unordered_set<int>& eids);

			void assignWeightsImpl(model* m, vector<float>& weights, const Config& config) override;
			unordered_set<int> m_cut_edges; //edges to be but using shortest path
		};
*/

	}
}
