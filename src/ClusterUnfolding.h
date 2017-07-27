/*
* ClusterUnfolding.h
*
*  Created on: Dec 01, 2015
*      Author: jmlien
*/

#pragma once

#include <string>
#include <vector>

#include "unfolder.h"
#include "Splitter.h"

namespace masc {
	namespace unfolding {

		using namespace std;

		struct FaceCluster
		{
			FaceCluster(){ seed = UINT_MAX;  total_score_error = 0; depth = 0; }
			list<uint> fids; //face ids
			uint seed; 
			uint depth;
			float total_score_error;
		};

		class ClusterUnfolding
		{
		public:

			ClusterUnfolding(Unfolder* unfolder);
			virtual ~ClusterUnfolding();


			bool setup(const string& filename);

			//////////////////////////////////////////
			// run n unfoldings and cluster
			//////////////////////////////////////////

			virtual void run();

			// repair the clusters from external labels
			void runRepair(const string& label_filename);

			//////////////////////////////////////////

			virtual void print(ostream& out) const;

			//////////////////////////////////////////
			// access
			//////////////////////////////////////////
			const vector<Unfolder *> getUnfolders() const {
			  vector<Unfolder *> output = this->m_successful_unfolded_meshes;
			  output.insert(output.end(), m_failed_unfolded_meshes.begin(), m_failed_unfolded_meshes.end());
			  return output;
			}

		protected:

			// will be called once after setup,
			// can be used to override parameters
			void virtual init();

			// setup the entire problem from stream
			virtual bool setup(istream& in);

			// setup the problem from tokens
			virtual bool setup(vector<string>& tokens);

			//run the heiristic k-times and collect k-unfoldings
			void simulate();

			// for each unfolding analyze wether a pair of faces can be unfolded 
			// considering (1) if the path connecting them overlap (2) if the path fit to the paper
			// this function will populate "m_score_matrix"
			void score();

			// estimate the score of faces on the model using the given splitter
			void score(const vector<float>& weight);

			//normalize the scoring matrix
			//find max score and compute new score = (1-score/max_score) 
			//so the best matching faces will be 0
			void normalize_scores(vector< vector<float> >& scores);

			//find max score and compute new score = (1-score/max_score)
			void normalize_scores(vector<float>& scores);

			// estimate the score of two faces on the model using the current unfolding
			// float score(uint f1, uint f2);

			//
			// create clusters
			void cluster(int cluster_size, vector<FaceCluster>& clusters);

			// cluster facets using the given seeds
			void cluster(unordered_set<uint>& seeds, vector<FaceCluster>& clusters);

			// compute a new seed from a cluster
			uint get_new_seed(const FaceCluster& cluster) const;

			// YH get neighbor triangle ids of a triangle
			vector<uint> get_neighbor_triangles(model * m, uint fid) const;
			vector<uint> get_neighbor_triangles(model * m, const triangle& t) const;

			// YH get neighbor (edge, triangle) ids of a triangle
			vector< pair<uint, uint> > get_neighbor_edges_triangles(model * mesh, uint fid) const;

			/// cluster the faces using spectral clustering
			void spectral_score_clustering(int cluster_size, vector< vector<float> >& scores, vector<FaceCluster>& clusters);


			/// cluster the faces using new spectral lib
			void spectral_clustering(int K, vector< vector<float> >& scores, vector<FaceCluster>& clusters);

			/// repair clusters so that each cluster contains a compact connected component
			void repair_clusters(model * m, vector<FaceCluster>& clusters, vector< vector<float> >& scores);


			/// get connected component of a give facet
			/// results are stored in cc
			void get_connected_component(model * m, uint fid, FaceCluster& cc) const;

			/// get connected components of a give cluster
			/// results are stored in ccs
			void get_connected_components(model * m, FaceCluster& cluster, list<FaceCluster>& ccs) const;

			// unfold a given cluster
			// return false, when unfolding failed
			// store unfolder if in the last iteration
			bool unfold_cluster(model * m, FaceCluster& cluster, bool last_iteration);

			/// build a submodel from the cluster
			model * build_mesh(model * m, FaceCluster & cluster);

			////check if there is collision/overlap along the path between the given fid and root
			//bool overlap_in_path_to_root(int fid, int root);

			/// find the first common ancestor of f1 and f2
			/// int first_common_ancestor(int f1, int f2);

			/// children points of the unfolded tree
			void process_unfolded_tree(vector< vector<uint> >& children, vector<uint>& overlap_faces) const;

			/// get vertices connected without passing through overlapping faces
			void getCC(uint root, vector<uint>& CC, const vector<vector<uint> >& children) const;

			/// save score to ppm image
			void save_score_to_image() const;
			
			/// save score to ppm image, using clustered order
			void save_score_to_image(vector<vector<uint> >& clusters) const;

			/// find N directions
			uint find_unfolding_directions(uint N, vector<Vector3d>& N_dirs);

			/// filename for serialization
			string get_serialization_filename() const;

			/// save score to file using filename from get_serialization_filename()
			void save_score_to_file() const;

			///load score from file using filename from get_serialization_filename()
			bool load_score_from_file();

			void dumpObj(const string& path, model * m) const;

		private:

			short count_same_connected_cluster(model * m, uint fid) const;

			void reassign_isolated_facet(model * m, uint fid, vector<FaceCluster>& clusters, const vector< vector<float> >& scores);

			void reassign_isolated_facets(model * m, list<uint>& fids, vector<FaceCluster>& clusters, const vector< vector<float> >& scores);

			void reassign_dangling_facet(model * m, uint fid, vector<FaceCluster>& clusters, const vector< vector<float> >& scores);

			void reassign_facet(model * m, uint fid, unordered_set<uint> adj_clusters, vector<FaceCluster>& clusters, const vector< vector<float> >& scores);

			double avg_dist_to_cluster(model * m, uint fid, FaceCluster& cluster, const vector< vector<float> >& scores);

			double cluster_area(model * m, const FaceCluster& cluster) const;

			void split_cluster(FaceCluster& cluster, int K, vector<FaceCluster>& new_clusters);

			Unfolder* m_unfolder;

			vector< vector<float> > m_weights; //edge weight of the mesh

			uint m_cluster_size; //number of clusters

			uint m_runs; //number of clusters

			float m_paper_width, m_paper_height; //paper width/height

			vector<std::string> m_heuristic_methods;       //name of the heritstic method, see main.h for detail
			vector<uint> m_m_heuristic_runs; //how many times should we run the heuristic method
			vector<CutHeuristic> m_heuristic_types;

			vector<Config> m_unfolder_configs; //use this to try final unfolding

			long m_start_time;

			vector< vector<float> > m_score_matrix; //this sparse matrix is indexed by face id. Each elemenet is <face_id, score>
			                                        //this is a similarity matrix
			vector<float> m_score_vector; // score for each face

			vector<FaceCluster> m_clusters;

			vector<Unfolder *> m_successful_unfolded_meshes;
			vector<Unfolder *> m_failed_unfolded_meshes;

			//
			enum CLUSTERING_METHOD { SPECTRAL, FLOYLD };
			CLUSTERING_METHOD m_clustering_method;

			//
			enum UNFOLD_DIRECTION { MODEL_NORMAL, RANDOM_DIRECTION };
			UNFOLD_DIRECTION m_unfold_dir;	

			unsigned char* m_visited;

			// only do clustering (without unfolding)
			// default is false
			bool m_clustering_only;

			//whether to diable repairing the cluster errors: isolated triangles
			// default is false
			bool m_disable_repair;
		};

	}
}
