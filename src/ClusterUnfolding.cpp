/*
* ClusterUnfolding.cpp
*
*  Created on: Dec 01, 2015
*      Author: jmlien
*/
#include "ClusterUnfolding.h"

#include <fstream>
#include <cassert>
#include <ctime>
#include <cfloat>
#include <algorithm>
#include <queue>    //YH
#include <sstream>
#include <string>
#include <list>
#include <unordered_set>
using namespace std;

#include <opencv2/core/core.hpp>

#include "spectralClustering.hpp"
#include "UnfoldingProblem.h"
#include "spatial_hash.h"
#include "objReader.h"

#include "libga/util/StreamHelper.h"
#include "util/UnfolderHelper.h"
#include "util/DataHelper.h"
#include "util/SVGWriter.h"

using namespace masc::util;
using namespace masc::unfolding::util;



#define threshold 1.0e-7
//#define num 30


namespace masc {
    namespace unfolding {

        const bool debug = false;
        const bool dump_cluster_obj=true;

        //defined at the bottom
        void normal_kmeans_clustering(int K, const vector<mathtool::Point3d> & all_dirs, vector<uint> & labels, vector<mathtool::Point3d> & centers);

        ClusterUnfolding::ClusterUnfolding(Unfolder* unfolder)
        {
            Config mycfg = unfolder->getConfig();
            mycfg.record_overlap = true;

            this->m_unfolder = new Unfolder(unfolder->getModel(), mycfg);
            assert(this->m_unfolder);
            m_cluster_size = 10;
            m_paper_width = m_paper_height = 100;

            m_clustering_method = SPECTRAL;  //SPECTRAL, FLOYLD
            m_unfold_dir = RANDOM_DIRECTION; //MODEL_NORMAL, RANDOM_DIRECTION
            m_runs = 0;

            m_visited = new unsigned char[unfolder->getModel()->t_size];
            assert(this->m_visited);

            m_start_time = 0;

            m_clustering_only = false; //default is to unfold the cluster

			m_disable_repair = false; //default is to repair the clustering error
        }

        ClusterUnfolding::~ClusterUnfolding()
        {
            delete this->m_unfolder;
            delete this->m_visited;
        }

        bool ClusterUnfolding::setup(const string& filename)
        {
            ifstream in;
            in.open(filename, ios_base::in);
            if (!in.good())
            {
                cerr << "!Error in ClusterUnfolding::setup. Can not open file = " << filename << endl;
                return false;
            }

            // invoke implementation
            auto flag = this->setup(in);

            in.close();

            return flag;
        }

        bool ClusterUnfolding::setup(istream& in)
        {
            const auto& org_config = this->m_unfolder->getConfig();

            while (true)
            {
                vector<string> tokens = ga::util::StreamHelper::tokenize(in);

                if (tokens.size() == 0) break;

                if (tokens[0] == "problem")
                {
                    if (!this->setup(tokens))
                    {
                        return false;
                    }

                    if (org_config.k > 0)
                    {
                        this->m_cluster_size = org_config.k;
                        cout << "!! cluster_size override to " << this->m_cluster_size;
                    }
                }
                else if (tokens[0] == "heuristic")
                {
                    for (int i = 0; i < tokens.size(); ++i)
                    {
                        const auto& token = tokens[i];

                        if (token == "method") {
                            this->m_heuristic_methods.push_back(tokens[++i]);
                        }
                        else if (token == "runs") {
                            this->m_m_heuristic_runs.push_back(stoi(tokens[++i]));
                        }
                    }
                }
                else if (tokens[0] == "unfold")
                {
                    for (int i = 0; i < tokens.size(); ++i)
                    {
                        const auto& token = tokens[i];

                        if (token == "method") {
                            const auto& name = tokens[++i];
                            Config config;
                            if (name == "GA" || name == "ga")
                            {
                                config.heuristic = CutHeuristic::GA;
                                config.ga_max_gen = org_config.ga_max_gen;
                                config.early_stop = org_config.early_stop;
                            }
                            m_unfolder_configs.push_back(config);
                        }
                        else if (token == "filename")
                        {
                            Config & config = m_unfolder_configs.back();
                            config.filename = config.ga_config_filename = tokens[++i];
                        }
                        else if (token == "seed")
                        {
                            Config & config = m_unfolder_configs.back();
                            config.seed = atoi(tokens[++i].c_str());
                        }
                    }
                }
                else if (tokens[0] == "paper")
                {
                    for (int i = 0; i < tokens.size(); ++i)
                    {
                        const auto& token = tokens[i];

                        if (token == "width") {
                            this->m_paper_width = stof(tokens[++i]);
                        }
                        else if (token == "height") {
                            this->m_paper_height = stof(tokens[++i]);
                        }
                    }
                }
                else
                {
                    return false;
                }
            }

            this->init();

            //srand48(time(nullptr));

            this->print(cout);


            return true;
        }

        bool ClusterUnfolding::setup(vector<string>& tokens)
        {
            //this->setSeed(time(nullptr));

            for (int i = 0; i<tokens.size(); ++i)
            {
                const auto& token = tokens[i];

                if (token == "size") {
                    this->m_cluster_size = stoi(tokens[++i]);
                }
                else if (token == "cluster")
                {
                    string name = tokens[++i];
                    if (name == "spectral")
                        m_clustering_method = SPECTRAL;
                }
                else if (token == "direction")
                {
                    string name = tokens[++i];
                    if (name == "random")
                        m_unfold_dir = RANDOM_DIRECTION;
                    else if (name == "model")
                        m_unfold_dir = MODEL_NORMAL;
                    //otherwise default is used (MODEL_NORMAL)
                }
				else if (token == "clustering-only")
					m_clustering_only = true;
				else if (token == "disable-repair")
					m_disable_repair = true;
            }
            return true;
        }

        // will be called once after setup,
        // can be used to override parameters
        void ClusterUnfolding::init()
        {
            for (string& method : m_heuristic_methods)
            {
                if (method[0] == 's')
                    m_heuristic_types.push_back(CutHeuristic::STEEPEST_EDGE);
                else if (method[0] == 'f')
                    m_heuristic_types.push_back(CutHeuristic::FLAT_TREE);
                else if (method[0] == 'p')
                    m_heuristic_types.push_back(CutHeuristic::MINIMUM_PERIMETER);
                else if (method[0] == 'r')
                    m_heuristic_types.push_back(CutHeuristic::RANDOM);
                else if (method[0] == 'b')
                    m_heuristic_types.push_back(CutHeuristic::BRUTE_FORCE);
                else //default
                    m_heuristic_types.push_back(CutHeuristic::STEEPEST_EDGE);
            }
        }

        void ClusterUnfolding::print(ostream& out) const
        {
            //nothing yet
        }


        void ClusterUnfolding::runRepair(const string& label_filename)
        {
          cout<<"- [ClusterUnfolding] Repairing mode..."<<endl;
          m_start_time = clock();
          const auto config = this->m_unfolder->getConfig();
          const auto m = this->m_unfolder->getModel();

          // pre allocate space
          m_score_matrix = vector<vector<float>>(m->t_size, vector<float>(m->t_size, 0));
          // load score
          this->load_score_from_file();

          // read the labels
          auto labels = masc::util::readList<int>(label_filename);

          cout<<"- [ClusterUnfolding] labels loaded, t_size = "<<labels.size()<<endl;

          assert(labels.size() == m->t_size);

          this->m_cluster_size = *std::max_element(labels.begin(), labels.end()) + 1;

          this->m_clusters = vector<FaceCluster>(m_cluster_size);

          for(int i=0;i<m->t_size;++i)
          {
            int cluster_id = labels[i];

            this->m_clusters[cluster_id].fids.push_back(i);
          }

          auto start = clock();
          this->repair_clusters(m, m_clusters, m_score_matrix);

          cout << "- [ClusterUnfolding] repairing time=" << (float)(clock() - start) / CLOCKS_PER_SEC << endl;

          // save repaired to file
          string repaired_label_filename = label_filename + ".repaired.seg";
          vector<int> cluster_ids;
          for(int i=0;i<m->t_size;++i)
            cluster_ids.push_back(m->tris[i].cluster_id);
          masc::util::writeList(repaired_label_filename, cluster_ids);

          cout << "- [ClusterUnfolding] repaired labels output to "<< repaired_label_filename << endl;
        }

        void ClusterUnfolding::run()
        {
          const auto config = this->m_unfolder->getConfig();

            m_start_time = clock();

            //run the heiristic k-times and collect k-unfoldings
            simulate();

            //
            //for each unfolding analyze wether a pair of faces can be unfolded
            //considering (1) if the path connecting them overlap (2) if the path fit to the paper
            //
            score();

            // only score, not do cluster
            if(config.score_only) return;

            //
            //the product of this step is an |F|X|F| matrix, each element is a measure of likelihood that
            //this pair of faces will be clustered
            //

            bool successful = false;

            auto mesh = this->m_unfolder->getModel();

            //create clusters
            cout << "- [ClusterUnfolding]  start clustering facets" << endl;
            if (this->m_clustering_method == FLOYLD)
                cluster(this->m_cluster_size, this->m_clusters);
            else
            {
                spectral_score_clustering(this->m_cluster_size, this->m_score_matrix, this->m_clusters);
				if (!m_disable_repair) repair_clusters(mesh, this->m_clusters, this->m_score_matrix);
            }

            int failed_size = 0;

            cout << "- [ClusterUnfolding]  There are " << m_clusters.size() << " clusters found" << endl;
#if 1
            vector<FaceCluster> todo, done;
            todo.swap(m_clusters);
            const int max_depth = 3;
            const int min_cluster_size = 200;
            const int split_size = 2;
            int iteration = 0;


            if (this->m_clustering_only == false) //skip unfolding if this is clustering only
            {

                while (todo.empty() == false && ++iteration <= config.max_iterations)
                {

                    cout << "- [ClusterUnfolding] iteration=" << iteration << " config.max_iterations=" << config.max_iterations << endl;
                    //unfold each cluster
                    uint cluster_id = 0;
                    vector<FaceCluster> failed;
                    for (auto& cluster : todo)
                    {
                        bool r = unfold_cluster(mesh, cluster, iteration == config.max_iterations);

                        //if there is a cluster failed to unfold, then increase the number of cluster by one
                        if (!r)
                        {
                            cout << "-------------------------------------------" << endl;
                            cout << "!!!!! Failed to unfold a subcomponent...." << endl;
                            cout << "-------------------------------------------" << endl;
                            failed.push_back(cluster);
                        }
                        else
                        {
                            cout << "-------------------------------------------" << endl;
                            cout << "- Great. This cluster is unfolded!! Move on" << endl;
                            cout << "-------------------------------------------" << endl;
                            done.push_back(cluster);
                        }
                    }

                    if (failed.empty())
                    {
                        //done
                        break;
                    }

                    todo.clear();

                    //handle failed clusters
                    cout << "[*** Failed size=" << failed.size() << " ***]" << endl;
                    for (auto& cluster : failed)
                    {
                        if (cluster.depth < max_depth && cluster.fids.size()>min_cluster_size) //max_depth and min_cluster_size is hard coded in this functon
                        {
                            //split the cluster and put then into todo
                            vector<FaceCluster> new_clusters;
                            split_cluster(cluster, split_size, new_clusters);
                            for (auto & nc : new_clusters)
                            {
                                nc.depth = cluster.depth + 1;
                                todo.push_back(nc);
                            }
                        }
                        else //report failure
                        {
                            auto subm = build_mesh(mesh, cluster);
                            stringstream ss;
                            ss << mesh->name << "_FAILED_unfolding_submodel_" << subm->t_size << "_" << subm->e_size << "_" << subm->v_size << ".obj";
                            this->dumpObj(ss.str(), subm);

							Config mycfg = m_unfolder->getConfig();
							mycfg.heuristic = CutHeuristic::STEEPEST_EDGE;
							mycfg.quite = true;
							mycfg.max_retries = 1;
							Unfolder * sub_unfolder = new Unfolder(subm, mycfg);
                            assert(sub_unfolder);
							sub_unfolder->buildUnfolding();
                            m_failed_unfolded_meshes.push_back(sub_unfolder);

							//subm->destroy();
							//delete subm;
							failed_size++;
                        }
                    }
                    //done handling failed clusters

                }//end of while
            }
            else{ //stuff for clustring only, put everything in succceful list
                Config mycfg = m_unfolder->getConfig();
                mycfg.heuristic = CutHeuristic::STEEPEST_EDGE;
                mycfg.quite = true;
                mycfg.max_retries = 1;
                for (auto& cluster : todo)
                {
                    auto subm = build_mesh(mesh, cluster);
                    Unfolder * sub_unfolder = new Unfolder(subm, mycfg);
                    assert(sub_unfolder);
                    //get a dummy unfolder
                    sub_unfolder->buildUnfolding();
                    m_successful_unfolded_meshes.push_back(sub_unfolder);
                    done.push_back(cluster);
                }
            }

            done.swap(this->m_clusters);
#endif

            if (failed_size == 0)
            {
                cout << "-------------------------------------------" << endl;
                cout << "- Great. All clusters are unfolded successfully!!" << endl;
                cout << "-------------------------------------------" << endl;
            }
            else
            {
                cout << "-------------------------------------------" << endl;
                cout << "- Failed to unfold " << failed_size<<" clusters..." << endl;
                cout << "-------------------------------------------" << endl;
            }

            cout << "- [ClusterUnfolding] this->m_successful_unfolded_meshes size=" << this->m_successful_unfolded_meshes.size() << endl;
            int cluster_id = 0;
            for (auto * unfolder : this->m_successful_unfolded_meshes)
            {
                unfolder->setClusterId(cluster_id);
                {
                    stringstream ss;
                    auto subm = unfolder->getModel();
                    ss << "./" << mesh->name << "_cluster_" << cluster_id << "_unfolding_submodel_" << subm->t_size << "_" << subm->e_size << "_" << subm->v_size << ".obj";
                    this->dumpObj(ss.str(), subm);
                }

                //{
                //	{
                //		stringstream ss;
                //		ss << "./" << mesh->name << "_cluster_" << cluster_id;
                //		unfolder->dumpSVG(ss.str() + ".svg"); //with everything
                //		unfolder->dumpSVG(ss.str() +"_cut.svg", 1); //1 boundary + creases
                //		unfolder->dumpSVG(ss.str() +"_extra.svg", 9); // 1 + extra cuts
                //		unfolder->dumpOri(ss.str() + ".ori"); // dump folding model
                //	}
                //}

                cluster_id++;
            }

            for(auto * unfolder : this->m_failed_unfolded_meshes)
            {
                unfolder->setClusterId(cluster_id++);
            }


            //-----------------------------------------------------------------
            /// record the cluster_id
            /// this may not be needed
            cluster_id = 0;
            for (auto& cluster : this->m_clusters)
            {
                for (auto& fid : cluster.fids)
                    mesh->tris[fid].cluster_id = cluster_id;
                cluster_id++;
            }
            //-----------------------------------------------------------------

        }// END ClusterUnfolding::run()


        void ClusterUnfolding::simulate()
        {
            int id = 0;
            mathtool::srand48(time(NULL));
            auto mesh = this->m_unfolder->getModel();

            for (auto type : this->m_heuristic_types)
            {
                //create a spliter
                auto splitter = Splitter::createSplitter(type);
                splitter->measure(mesh);
                uint runs = this->m_m_heuristic_runs[id++];

                //
                if (type == CutHeuristic::STEEPEST_EDGE || type == CutHeuristic::FLAT_TREE) //YH default
                {
                    //create random directions (this block of code will be modified so the directions are collected from the model)
                    vector<Vector3d> rand_dirs;

                    if (this->m_unfold_dir == RANDOM_DIRECTION)
                    {

                        for (uint i = 0; i < runs; i++)
                        {
                            Vector3d rv(mathtool::drand48(), mathtool::drand48(), mathtool::drand48());
                            rv = rv.normalize();
                            if (mathtool::drand48()>0.5) rv[0] = -rv[0];
                            if (mathtool::drand48() > 0.5) rv[1] = -rv[1];
                            if (mathtool::drand48() > 0.5) rv[2] = -rv[2];
                            rand_dirs.push_back(rv);
                        }
                    }
                    else //use normal directions from the model
                    {
                        runs = find_unfolding_directions(runs, rand_dirs);
                    }

                    //ask splitter to generated weight
                    Config mycfg=m_unfolder->getConfig();
                    mycfg.use_user_vector = true;
                    mycfg.quite = true;
                    mycfg.early_stop=true;
                    for (int r = 0; r < runs; r++)
                    {
                        mycfg.user_vector = rand_dirs[r];
                        this->m_weights.push_back(splitter->assignWeights(mesh, mycfg));
                    }
                }
                else if (type == CutHeuristic::RANDOM)
                {
                    Config mycfg = m_unfolder->getConfig();
                    mycfg.quite = true;
                    for (int r = 0; r < runs; r++)
                        this->m_weights.push_back(splitter->assignWeights(mesh, mycfg));
                }
                else //these are determinisitic methods
                {
                    Config mycfg = m_unfolder->getConfig();
                    mycfg.quite = true;
                    this->m_weights.push_back(splitter->assignWeights(mesh, mycfg));
                }

                delete splitter;
            }//end for type
        }

        void ClusterUnfolding::score()
        {
            //fill in m_score_matrix
            const auto config = this->m_unfolder->getConfig();
            uint size = this->m_unfolder->getModel()->t_size;
            if(config.scalar_score) {
              this->m_score_vector = vector<float>(size, 0);
            } else {
              this->m_score_matrix = vector< vector<float> >(size, vector<float>(size, 0));
            }


            if (!config.force_mode && load_score_from_file())
            {
                cout << "- [ClusterUnfolding] loaded " << this->m_runs<<" from file : " << get_serialization_filename() << endl;
                int remain = this->m_weights.size() - this->m_runs;
                if(remain>0) cout << "- still need " << remain << " runs";
                else remain = 0;
                m_weights.erase(m_weights.begin() + remain, m_weights.end());
                //auto tmp_score_matrix = this->m_score_matrix;
                //this->m_score_matrix = vector< vector<float> >(size, vector<float>(size, 0));

                //for (auto& weight : this->m_weights)
                //{
                //	score(weight);
                //}

                ////normalize the scoring matrix
                ////find max score and compute new score = (1-score/max_score)
                ////so the best matching faces will be 0
                //normalize_scores();

                ////compare....
                //for (int i = 0; i < this->m_unfolder->getModel()->t_size; i++)
                //{
                //	for (int j = 0; j < this->m_unfolder->getModel()->t_size; j++)
                //	{
                //		auto diff = fabs(tmp_score_matrix[i][j] - this->m_score_matrix[i][j]);
                //		if (diff>1e-10)
                //			cout << "!!!!!!!!!!!!!!!!!!!!!! help! [" << i << "][" << j << "] diff=" << diff << " loaded data=" << tmp_score_matrix[i][j]<<endl;
                //	}
                //}
            }
            else
                this->m_runs = 0; //nothing is loaded

            if (m_weights.empty()==false)
            {
                int total_size = this->m_weights.size();
                int counter = 0;
                for (auto& weight : this->m_weights)
                {
                    cout << "- [ClusterUnfolding] scoring [" << counter++ << "/" << total_size << "]" <<"\r"<<flush;
                    score(weight);
                }
                this->m_runs += total_size;
                cout << "\n- [ClusterUnfolding] scoring done" << endl;

                //save score to file
                save_score_to_file();
            }

            if(config.scalar_score)
            {
              normalize_scores(this->m_score_vector);
              save_score_to_file();

              for(int i=0;i<size;++i)
                this->m_unfolder->getModel()->tris[i].score = this->m_score_vector[i];
            }
            else
            {
                //normalize the scoring matrix
                //find max score and compute new score = (1-score/max_score)
                //so the best matching faces will be 0
                normalize_scores(this->m_score_matrix);
            }

            //for debug only
            if (debug)
                save_score_to_image();
        }

        /*
        //estimate the score of faces on the model using the given splitter
        void ClusterUnfolding::score(const vector<float>& weight) {
            //create a new unfolder...
            {
                Config mycfg = m_unfolder->getConfig();
                mycfg.quite = true;
                Unfolder * new_unfolder = new Unfolder(m_unfolder->getModel(), mycfg);
                assert(new_unfolder);
                delete m_unfolder;
                m_unfolder = new_unfolder;
            }

            auto start = clock();

            int overlap_count = this->m_unfolder->buildFromWeights(weight);
            this->m_unfolder->rebuildModel();

            if (debug) {

                static int count = 0;
                count++;
                stringstream ss;
                ss << this->m_unfolder->getModel()->name << "_unfolding_" << count
                    << ".svg";
                this->m_unfolder->dumpSVG(ss.str(), 1);

                cout << "[" << count << "] overlap_count = " << overlap_count << endl;
            }

            if (overlap_count == 0) {

                for (auto& f_scores : this->m_score_matrix) {
                    for (auto& sc : f_scores) {
                        sc += 1;
                    }			//end for f_scores
                }
                return; //a valid unfolding is found!?
            }

            model * m = m_unfolder->getModel();

            //this->m_unfolder->rebuildModel();
            vector<vector<uint> > children;
            vector<uint> overlap_faces;
            process_unfolded_tree(children, overlap_faces);

            cout << "unfold time=" << (float)(clock() - start) / CLOCKS_PER_SEC << endl;
            cout << "overlap_faces size=" << overlap_faces.size() << endl;

            //get subtree roots
            //unordered_set<uint> roots;
            //{
            //	roots.insert(this->m_unfolder->getBaseFaceID());
            //	for (uint of : overlap_faces)
            //	{
            //		list<uint> of_kids = children[of];
            //		roots.insert(of);
            //		for (uint kid : of_kids)
            //		{
            //			roots.insert(kid);
            //		}
            //	}
            //}

            //get CC (connected components)
            //uint root_size = roots.size();
            //for (uint root : roots)

            auto& collisions = m_unfolder->getOverlppingFacePairs();

            start = clock();
            uint root_size = m->t_size;
            vector<uint> CC;
            for (uint root = 0; root < root_size; root++) {
                CC.clear();
                getCC(root, CC, children);

                //update score matrix
                uint cc_size = CC.size();
                for (int i = 0; i < cc_size; i++) {
                    uint id = CC[i];
                    if (id == root)
                        continue;
                    this->m_score_matrix[root][id] += (1.0);
                    this->m_score_matrix[id][root] += (1.0);
                }
            }
            cout << "score time=" << (float)(clock() - start) / CLOCKS_PER_SEC << endl;
        }
        */

        //estimate the score of faces on the model using the given splitter
        void ClusterUnfolding::score(const vector<float>& weight)
        {
            model * m = m_unfolder->getModel();
            const auto config = m_unfolder->getConfig();

            //create a new unfolder...
            {
                Config mycfg = m_unfolder->getConfig();
                mycfg.quite = true;
                Unfolder * new_unfolder = new Unfolder(m, mycfg);
                assert(new_unfolder);
                delete m_unfolder;
                m_unfolder = new_unfolder;
            }

            auto start = clock();

            int overlap_count = this->m_unfolder->buildFromWeights(weight);
            this->m_unfolder->rebuildModel();

            if (debug)
            {
                static int count = 0;
                count++;
                stringstream ss;
                ss << "-  [ClusterUnfolding] "<<m->name << "_unfolding_" << count << ".svg";



                m_unfolder->dumpSVG(ss.str(), ExportSVGType::NORMAL);

                cout << "- [ClusterUnfolding] [" << count << "] overlap_count = " << overlap_count << endl;
            }

            UnfolderHelper unfolder_helper(m_unfolder);


            if (overlap_count == 0)
            {
                // scalar score
                if (config.scalar_score)
                {
                    // add entire surface area to all faces
                    for (int i = 0; i < m->t_size; ++i)
                        this->m_score_vector[i] += m->surface_area;
                }
                else
                {
                    if (config.weighted_dist)
                    {
                        auto dist = unfolder_helper.computeGeoDist();
                        for (int i = 0; i < m->t_size; ++i)
                        for (int j = i + 1; j < m->t_size; ++j)
                        {
                            this->m_score_matrix[i][j] += 1.0 / dist[i][j];
                            this->m_score_matrix[j][i] += 1.0 / dist[i][j];
                        }
                    }
                    else {
                        for (auto& f_scores : this->m_score_matrix)
                        {
                            for (auto& sc : f_scores)
                            {
                                sc += 1;
                            }//end for f_scores
                        }
                    }
                } // end if scalar_score

                return; //a valid unfolding is found!?
            }



            //this->m_unfolder->rebuildModel();
            vector< vector<uint> > children;
            vector<uint> overlap_faces;
            process_unfolded_tree(children, overlap_faces);

      			if (!config.quite)
      			{
      				cout << "- [ClusterUnfolding] unfold time=" << (float)(clock() - start) / CLOCKS_PER_SEC << endl;
      				cout << "- [ClusterUnfolding] overlap_faces size=" << overlap_faces.size() << endl;
      			}

            //get subtree roots
            //unordered_set<uint> roots;
            //{
            //	roots.insert(this->m_unfolder->getBaseFaceID());
            //	for (uint of : overlap_faces)
            //	{
            //		list<uint> of_kids = children[of];
            //		roots.insert(of);
            //		for (uint kid : of_kids)
            //		{
            //			roots.insert(kid);
            //		}
            //	}
            //}

            //get CC (connected components)
            //uint root_size = roots.size();
            //for (uint root : roots)


            auto& collisions = m_unfolder->getOverlppingFacePairs();

            start = clock();
            uint root_size = m->t_size;

            UnfolderHelper helper(this->m_unfolder);

            unordered_set<uint> visited_fids;

            for (uint root = 0; root<root_size; root++)
            {
              if(config.scalar_score && visited_fids.count(root)) {
                // already found a cc that contains current face...
                continue;
              }

                vector<uint> CC;
                getCC(root, CC, children);

                //update score matrix
                uint cc_size = CC.size();

                if(config.scalar_score) {
                  float unfolding_area = 0.0f;
                  // scalar score
                  for(auto fid : CC)
                  {
                    unfolding_area += (m->tris[fid].area);
                    visited_fids.insert(fid);
                  }
                  for(auto fid : CC)
                    this->m_score_vector[fid] += (unfolding_area - m->tris[fid].area);
                }
                else {
                  // pairwise score
          if(config.weighted_dist)
          {
            auto dist = helper.computeGeoDist(root, CC);
            for(int i=0;i<cc_size;++i) {
              uint id = CC[i];
              if(id == root) continue;
              //if(dist[i]==0.0) {
              //  cout<<"dist == 0 root = "<<root<<" id = "<<id<<endl;
              //}
              this->m_score_matrix[root][id] += (1.0 / dist[i]);
              this->m_score_matrix[id][root] += (1.0 / dist[i]);
            }
          }
          else {
            for (int i = 0; i < cc_size; i++)
            {
              uint id = CC[i];
              if (id == root) continue;
              this->m_score_matrix[root][id] += (1.0);
              this->m_score_matrix[id][root] += (1.0);
            }
          }
                }
            }
            //cout << "score time=" << (float)(clock() - start) / CLOCKS_PER_SEC << "\r";
        }

        void ClusterUnfolding::normalize_scores(vector<float>& scores)
        {
          const auto max_score = *std::max_element(scores.begin(), scores.end());
          for(auto& score : scores)
            score /= max_score;
          const auto new_max_score = *std::max_element(scores.begin(), scores.end());
          cout<<"-  [ClusterUnfolding] max_score = "<<max_score<<"new max_score"<<new_max_score<<endl;
        }

        //nomrlize scores in the matrix
        void ClusterUnfolding::normalize_scores(vector< vector<float> >& scores)
        {
            //normalize the scoring matrix
            //find max score and compute new score = (1-score/max_score)
            //so the best matching faces will be 0
            float max_score = 0;
            float min_score = FLT_MAX;
            uint min_row = 0;
            uint size = scores.size();
            for (int i = 0; i<size; i++)
            {
                float row_score = 0;
                for (int j = 0; j<size; j++)
                {
                    auto sc = scores[i][j];
                    if (sc > max_score) max_score = sc;
                    row_score += sc;
                }

                if (row_score < min_score){
                    min_score = row_score;
                    min_row = i;
                }

            }////end for f_scores

            cout << "- [ClusterUnfolding] max score = " << max_score << " min_row=" << min_row << " min_score = " << min_score << endl;


            for (int i = 0; i<size; i++)
            {
                for (int j = 0; j<size; j++)
                {
                    auto& sc = scores[i][j];
                    //if(i!=j)
                    sc = (1 - sc / max_score);
                    //if(i==j) the score remains 0
                }
            }////end for f_scores
        }

        //
        // estimate the score of two faces on the model using the current unfolding
        // 0 means that f1 and f2 are connected on the unfolding without overlapping
        // 1 means that f1 and f2 are connected via overlapping faces
        //
        //float ClusterUnfolding::score(uint f1, uint f2)
        //{
        //	int common_ancestor=first_common_ancestor(f1, f2);
        //	if (overlap_in_path_to_root(f1, common_ancestor)) return 1;
        //	if (overlap_in_path_to_root(f2, common_ancestor)) return 1;
        //	return 0;
        //}

        void ClusterUnfolding::cluster(int cluster_size, vector<FaceCluster>& clusters)	//YH this is my part
        {
            ////print cout
            //cout << "===========================================" << endl;
            //cout << "============= CLUSTERUNFOLDING ============" << endl;
            //cout << "===========================================" << endl;

            //mathtool::srand48(time(NULL));


            ////get a mesh model
            //auto mesh = m_unfolder->getModel();

            ////print cout
            //cout << "  Cluster size: " << cluster_size << endl;
            //cout << "  Face size: " << mesh->t_size << endl;
            //
            ////set seed face id; at first, set ids randomly
            //unordered_set<uint> seed_fids;

            ////set face id randomly
            //if (cluster_size > mesh->t_size)
            //{
            //	cout << "Error! cluster size is too big" << endl;
            //	exit(-1);
            //}
            ////else
            ////{
            ////	while (seed_fids.size()<cluster_size)
            ////	{
            ////		uint fid = drand48()*mesh->t_size; //create random face id
            ////		seed_fids.insert(fid);
            ////	}
            ////}

            ////spectral_score_clustering(this->m_cluster_size, this->m_clusters);
            ////for (auto& cluster : clusters)
            ////{
            ////	auto seed = this->get_new_seed(cluster);
            ////	seed_fids.insert(seed);
            ////}

            //cout << "  Random face id: ";
            //for (auto a : seed_fids)	cout << a << "--";
            //cout << endl;

            //bool flag = true; //finish iteration
            ////int num_iterator = 0;
            //double previous_total_errors = DBL_MAX;

            //while (flag) //given seed fids, calculate their scores
            //{
            //	clusters.clear();
            //	cluster(seed_fids, clusters);

            //	double total_errors=0;
            //	int new_seed_cout = 0;
            //	{//for each cluster, find another seed face ids which has the lowest score in the cluster
            //		seed_fids.clear();
            //		for (auto& cluster : clusters)
            //		{
            //			cout << "cluster size=" << cluster.fids.size() << endl;
            //			total_errors += cluster.total_score_error;
            //			auto seed = this->get_new_seed(cluster);
            //			seed_fids.insert(seed);
            //			if (seed != cluster.seed) new_seed_cout++;
            //		}
            //	}

            //	cout << "total_errors=" << total_errors << endl;
            //	{//check if we need to iterate
            //		if (new_seed_cout==0)
            //		{
            //			cout << "terminate: all seeds remain the same" << endl;
            //			flag = false;
            //		}

            //		if (fabs(total_errors - previous_total_errors) < threshold)
            //		{
            //			cout << "terminate: converged" << endl;
            //			flag = false;
            //		}
            //	}

            //	previous_total_errors = total_errors;

            //	break;

            //}

        }

        void ClusterUnfolding::cluster(unordered_set<uint>& seeds, vector<FaceCluster>& clusters)	//YH this is my part
        {
            //struct cluster_assignment
            //{
            //	cluster_assignment(){ seed = UINT_MAX;  score = FLT_MAX; }
            //	uint fid; //facet id
            //	uint seed; //seed face id
            //	float score; //score of the face fid assigned to cluster cid

            //	bool operator<(const cluster_assignment& other) const {
            //		return score > other.score;
            //	}
            //};


            //struct face_assignment
            //{
            //	face_assignment(){ seed = UINT_MAX;  score = FLT_MAX; }
            //	uint seed; //seed face id
            //	float score; //score of the face fid assigned to cluster cid
            //};
            //vector<face_assignment> f_ass(m_unfolder->getModel()->t_size, face_assignment());


            ////priority queue along the scores
            ////put all face data in q later
            //vector<cluster_assignment> q;

            //{//push the seed face ids to q

            //	for (auto id : seeds)
            //	{
            //		cluster_assignment tmp_ca;

            //		tmp_ca.fid = id;
            //		tmp_ca.seed = id;
            //		tmp_ca.score = this->m_score_matrix[id][id];
            //		q.push_back(tmp_ca);
            //	}

            //	//sort q along the score
            //	make_heap(q.begin(), q.end());
            //	//sort_heap(q.begin(), q.end());
            //}

            //cout << "  -- initial q: " << endl;
            //for (auto a : q) cout << "  fid: " << a.fid << " seed: " << a.seed << " score: " << a.score << endl;

            //while (!q.empty())
            //{
            //
            //	auto first_q = q.front();
            //	pop_heap(q.begin(), q.end());
            //	q.pop_back();

            //	if (f_ass[first_q.fid].score <= first_q.score) continue; //nothing to look at here

            //	//smaller score is found
            //	f_ass[first_q.fid].score = first_q.score;
            //	f_ass[first_q.fid].seed = first_q.seed;

            //	//start to propagate
            //	vector<uint> nei_fids = get_neighbor_triangles(m, first_q.fid); //get neighbors of first_q.fid

            //	for (auto first_nei_fids : nei_fids)
            //	{
            //		//float score = first_q.score + this->m_score_matrix[first_q.seed][first_nei_fids];
            //		float score = this->m_score_matrix[first_q.seed][first_nei_fids];
            //		//this->m_score_matrix[first_q.fid][first_nei_fids] +

            //		if (score < f_ass[first_nei_fids].score) //find a better score
            //		{
            //			cluster_assignment tmp_ac;

            //			tmp_ac.fid = first_nei_fids;
            //			tmp_ac.score = score;
            //			tmp_ac.seed = first_q.seed;

            //			q.push_back(tmp_ac);
            //			push_heap(q.begin(), q.end());
            //		}
            //	}//end for first_nei_fids

            //}//end while loop

            ////convert vector<face_assignment> f_ass to vector<FaceCluster>& clusters
            //map<uint, FaceCluster> seed2clusters;
            //uint tsize = m_unfolder->getModel()->t_size;
            //for (uint fid = 0; fid < tsize;fid++)
            //{
            //	auto & s2c = seed2clusters[f_ass[fid].seed];
            //	s2c.fids.push_back(fid);
            //	s2c.total_score_error += f_ass[fid].score;
            //}

            ////store
            //for (auto & s2c : seed2clusters)
            //{
            //	s2c.second.seed = s2c.first;
            //	clusters.push_back(s2c.second);
            //}

            ////done
        }

        inline double area(model * m, triangle& tri)
        {
            auto& v0 = m->vertices[tri.v[0]];
            auto& v1 = m->vertices[tri.v[1]];
            auto& v2 = m->vertices[tri.v[2]];
            Vector3d vec1 = v1.p - v0.p;
            Vector3d vec2 = v2.p - v0.p;
            return (vec1%vec2).norm() / 2;
        }

        uint ClusterUnfolding::get_new_seed(const FaceCluster& cluster) const	//YH this is my part
        {
            double min_dist = DBL_MAX;
            uint best_fid = cluster.seed;

            for (auto fid1 : cluster.fids)
            {
                double total=0;
                for (auto fid2 : cluster.fids)
                {
                    if (fid1 == fid2) continue;
                    total += this->m_score_matrix[fid1][fid2] * area(this->m_unfolder->getModel(), this->m_unfolder->getModel()->tris[fid2]);
                }

                if (total < min_dist)
                {
                    min_dist = total;
                    best_fid = fid1;
                }
            }

            return best_fid;
        }

        model * ClusterUnfolding::build_mesh(model * m, FaceCluster & cluster)
        {
            //check if there is non-overlapping unfolding for this cluster
            //
            //model * m = m_unfolder->getModel();
            objModel data;

            //collect all vertices
            unordered_map<uint, uint> vids;
            for (auto fid : cluster.fids)
            {
                for (short d = 0; d<3; d++) vids[m->tris[fid].v[d]] = 0;
            }

            //get a list of points
            uint new_vid = 0;
            for (auto & vid : vids)
            {
                Vpt pt;
                auto & pos = m->vertices[vid.first].p;
                pt.x = pos[0];
                pt.y = pos[1];
                pt.z = pos[2];
                data.pts.push_back(pt);
                vid.second = new_vid++;
            }

            //get a list of faces
            for (auto fid : cluster.fids)
            {
                polygon poly;
                for (short d = 0; d < 3; d++){
                    poly.pts.push_back(vids[m->tris[fid].v[d]]);
                }
                data.polys.push_back(poly);
            }

            data.compute_v_normal();


            //build mesh
            model * subm = new model();
            if (subm->build(data) == false)
            {
                cerr << "! Error: Failed to build a model from face cluster" << endl;
                return NULL;
            }

            //register
            auto fid_it = cluster.fids.begin();
            for (int fid = 0; fid < subm->t_size; fid++)
            {
                subm->tris[fid].source_fid = *fid_it;
                ++fid_it;
            }

            return subm;
        }

        void ClusterUnfolding::split_cluster(FaceCluster& cluster, int K, vector<FaceCluster>& new_clusters)
        {
            int fsize = cluster.fids.size();
            vector<uint> fids(cluster.fids.begin(), cluster.fids.end());
            model * subm = build_mesh(this->m_unfolder->getModel(), cluster);

            vector< vector<float> > sub_scores(fsize, vector<float>(fsize, 0));
            for (int x = 0; x < fsize;x++)
            {
                for (int y = 0; y < fsize; y++)
                {
                    sub_scores[x][y] = this->m_score_matrix[fids[x]][fids[y]];
                }
            }

            //normalize_scores(sub_scores);

            spectral_score_clustering(K, sub_scores, new_clusters);
			if (!m_disable_repair) repair_clusters(subm, new_clusters, sub_scores);

            for (auto& cluster : new_clusters)
            {
                list<uint> new_fids;
                for (uint fid : cluster.fids)
                {
                    new_fids.push_back(fids[fid]);
                }
                cluster.fids.swap(new_fids);
            }

            subm->destroy();
            delete subm;
        }

        bool ClusterUnfolding::unfold_cluster(model * m, FaceCluster& cluster, bool last_iteration)
        {
            cout << "- [ClusterUnfolding] last_iteration=" << last_iteration << endl;
            //return false;

            cout << "- [ClusterUnfolding] --------> start to unfold a cluster of size=" << cluster.fids.size() << " depth=" << cluster .depth << endl;

            bool successfull = false;

            //build mesh
            model * subm = build_mesh(m, cluster);
            if (subm == NULL)
            {
                cerr << "! [ClusterUnfolding] Error: Failed to build a model from face cluster" << endl;
                return false;
            }

            //now try to unfold this model
            uint heuristic_id = 0;

            for (auto & cfg : this->m_unfolder_configs)
            {
                // cfg.record_overlap = false;
                cfg.quite = true;
                cfg.early_stop=true; //stop when find a net

                Unfolder * sub_unfolder = new Unfolder(subm, cfg);

                if (cfg.heuristic == CutHeuristic::GA)
                {

                    // using ga to find the best unfolding
                    UnfoldingProblem problem(sub_unfolder);

                    if (problem.setup(cfg.ga_config_filename))
                    {
                        // don't forget to set seed since it use another random generator instead of drand48 function
                        problem.setSeed(cfg.seed);
                        problem.run();
                        if (problem.isGoalAchieved()) {
                          sub_unfolder->rebuildModel();
                          m_successful_unfolded_meshes.push_back(sub_unfolder);
                          return true;
                        }
                        else if(last_iteration) {
                          sub_unfolder->rebuildModel();
                          m_failed_unfolded_meshes.push_back(sub_unfolder);
                          return false;
                        }
                    }
                }
                else{
                    cerr << "! [ClusterUnfolding] Error: unknown unfolder" << endl;
                    exit(1);
                }

                delete sub_unfolder;
            }

            delete subm;
            return false; //failed...
        }

        // YH get neighbor triangle id
        vector<uint> ClusterUnfolding::get_neighbor_triangles(model * mesh, uint fid) const
        {
            vector<uint> nei_fids; //neighbor face ids
            triangle& t = mesh->tris[fid];

            for (int i = 0; i < 3; i++){
                edge& e = mesh->edges[t.e[i]];
                //vector<uint> f = e.fid;
                uint nei_fid = (fid == e.fid[0]) ? e.fid[1] : e.fid[0];

                nei_fids.push_back(nei_fid);
            }

            return nei_fids;
        }

        // YH get neighbor triangle id
        vector<uint> ClusterUnfolding::get_neighbor_triangles(model * mesh, const triangle& t) const
        {
            vector<uint> nei_fids; //neighbor face ids

            for (int i = 0; i < 3; i++){
                edge& e = mesh->edges[t.e[i]];
                //vector<uint>& f = e.fid;
                uint nei_fid = (&t == &(mesh->tris[e.fid[0]])) ? e.fid[1] : e.fid[0];

                nei_fids.push_back(nei_fid);
            }

            return nei_fids;
        }

        /// children points of the unfolded tree
        void ClusterUnfolding::process_unfolded_tree(vector< vector<uint> >& children, vector<uint>& overlap_faces) const
        {
            //cout << "-----------------" << endl;

            const model * m = m_unfolder->getModel();
            children.resize(m->t_size);
            for (uint i = 0; i < m->t_size; i++)
            {
                auto pid = m_unfolder->getParentFaceId(i);
                if (pid == -1) continue; //i is the root
                children[pid].push_back(i);

                if (m->tris[i].overlapped)
                {
                    overlap_faces.push_back(i);
                    //cout << "overlap fid=" << i << endl;
                }
            }
        }//end determine_children


        //get a list of non-overlapping faces rooted at "root"
        void ClusterUnfolding::getCC(uint root, vector<uint>& CC, const vector<vector<uint> >& children) const
        {

            const model * m = m_unfolder->getModel();
            const auto& collisions = m_unfolder->getOverlppingFacePairs();
            deque<uint> open;

            //record what has been included in this CC
            //  bool * visited = new bool[m->t_size];
            // memset(visited, false, sizeof(bool) * m->t_size);

            memset(m_visited, 0, sizeof(unsigned char)* m->t_size);

            //add root
            m_visited[root] = true;
            open.push_back(root);

            //loop until all faces in open are handled
            while (open.empty() == false) {
                uint fid = open.front();
                open.pop_front();
                CC.push_back(fid);

                vector<uint> next = children[fid];
                uint pid = m_unfolder->getParentFaceId(fid);
                if (pid != -1)
                    next.push_back(pid);

                for (const auto kid : next) {
                    if (m_visited[kid])
                        continue;

                    //a free face, simply add it to CC
                    if (m->tris[kid].overlapped == false) {
                        open.push_back(kid);
                        m_visited[kid] = true;
                    }
                    else //now that the kid is conflicting with another face.
                        // we need to make sure that it does not conflict with the faces in CC
                    {
                        const auto& conflicts = collisions[kid];
                        bool found_conflict = false;
                        for (const auto conflict : conflicts) {
                            if (m_visited[conflict]) //the conflict is already in CC
                            {
                                found_conflict = true;
                                break;
                            }
                        } //end for conflict

                        if (found_conflict == false) //no conflict found, good to go
                        {
                            open.push_back(kid);
                            m_visited[kid] = true;
                        }
                    }
                } //end for kid

            } //end while

            //  delete[] visited;
        }




        ////get a list of non-overlapping faces rooted at "root"
        //void ClusterUnfolding::getCC(uint root, vector<uint>& CC, const vector< list<uint> >& children) const
        //{

        //	const model * m = m_unfolder->getModel();
        //	auto& collisions = m_unfolder->getOverlppingFacePairs();
        //	list<uint> open;

        //	//record what has been included in this CC
        //	bool * visited = new bool[m->t_size];
        //	memset(visited, false, sizeof(bool)*m->t_size);

        //	//add root
        //	visited[root] = true;
        //	open.push_back(root);

        //	//loop until all faces in open are handled
        //	while (open.empty() == false)
        //	{
        //		uint fid = open.front();
        //		open.pop_front();
        //		CC.push_back(fid);

        //		list<uint> next = children[fid];
        //		uint pid = m_unfolder->getParentFaceId(fid);
        //		if (pid != -1) next.push_back(pid);

        //		for (auto kid : next)
        //		{
        //			if (visited[kid]) continue;

        //			//a free face, simply add it to CC
        //			if (m->tris[kid].overlapped == false)
        //			{
        //				open.push_back(kid);
        //				visited[kid] = true;
        //			}
        //			else //now that the kid is conflicting with another face.
        //				// we need to make sure that it does not conflict with the faces in CC
        //			{
        //				auto& conflicts = collisions[kid];
        //				bool found_conflict = false;
        //				for (auto conflict : conflicts)
        //				{
        //					if (visited[conflict]) //the conflict is already in CC
        //					{
        //						found_conflict = true;
        //						break;
        //					}
        //				}//end for conflict

        //				if (found_conflict == false) //no conflict found, good to go
        //				{
        //					open.push_back(kid);
        //					visited[kid] = true;
        //				}
        //			}
        //		}//end for kid

        //	}//end while

        //	delete[] visited;
        //}

        /// save score to ppm image
        void ClusterUnfolding::save_score_to_image() const
        {
            string filename = "unfolding_score_matrix.ppm";
            FILE *f = fopen(filename.c_str(), "w");         // Write image to PPM file.

            uint size = m_unfolder->getModel()->t_size;
            fprintf(f, "P3\n%d %d\n%d\n", size, size, 255);

            for (int i = 0; i < size; i++)
            {
                for (int j = 0; j < size; j++)
                {
                    //higher intensity means higher correlation
                    int intensity = (int)((1 - m_score_matrix[i][j]) * 255);
                    fprintf(f, "%d %d %d ", intensity, 0, 0);
                }

            }
            fclose(f);
        }

        /// save score to ppm image
        void ClusterUnfolding::save_score_to_image(vector<vector<uint> >& clusters) const
        {
            string filename = "clustered_unfolding_score_matrix.ppm";
            FILE *f = fopen(filename.c_str(), "w");         // Write image to PPM file.

            uint size = m_unfolder->getModel()->t_size;
            fprintf(f, "P3\n%d %d\n%d\n", size, size, 255);

            vector<uint> new_order;
            vector<uint> idx(size, 0);

            int cid = 0;//cluster id
            vector<Point3d> cluster_color;
            for (auto& cluster : clusters)
            {
                for (auto id : cluster)
                {
                    new_order.push_back(id);
                    idx[id] = cid;
                }

                switch (cid % 6)
                {
                case 0:cluster_color.push_back(Point3d(255, 0, 0)); break;
                case 1:cluster_color.push_back(Point3d(255, 0, 255)); break;
                case 2:cluster_color.push_back(Point3d(0, 255, 0)); break;
                case 3:cluster_color.push_back(Point3d(0, 0, 255)); break;
                case 4:cluster_color.push_back(Point3d(255, 255, 0)); break;
                case 5:cluster_color.push_back(Point3d(0, 255, 255)); break;
                }
                cid++;
            }

            for (int i = 0; i < size; i++)
            {
                uint id1 = new_order[i];
                for (int j = 0; j < size; j++)
                {
                    uint id2 = new_order[j];

                    if (idx[id1] == idx[id2])
                    {
                        const auto& color = cluster_color[idx[id1]];
                        float intensity = (1 - m_score_matrix[id1][id2]);
                        fprintf(f, "%d %d %d ", (int)(intensity*color[0]), (int)(intensity*color[1]), (int)(intensity*color[2]));
                        //int intensity = (int)((1 - m_score_matrix[id1][id2]) * 255);
                        //fprintf(f, "%d %d %d ", intensity, 0, 0);
                    }
                    else
                    {
                        //higher intensity means higher correlation
                        int intensity = (int)((1 - m_score_matrix[id1][id2]) * 255);
                        fprintf(f, "%d %d %d ", intensity, intensity, intensity);
                    }
                }

            }
            fclose(f);
        }

        template<typename V> bool isfinite(V& vec)
        {
            if (std::isfinite(vec[0]) == false) return false;
            if (std::isfinite(vec[1]) == false) return false;
            if (std::isfinite(vec[2]) == false) return false;
            return true;
        }
        /// find N directions
        uint ClusterUnfolding::find_unfolding_directions(uint N, vector<Vector3d>& N_dirs)
        {
            if (N == 0) return 0; //nothing to do

            const model * m = m_unfolder->getModel();

            vector<mathtool::Point3d> all_dirs;
            {
                {
                    spatial_hash hash(mathtool::Point3d(-1, -1, -1), 2, 2, 2);

                    for (uint t = 0; t < m->t_size; t++) {
                        //this is necessary bc some faces have zero areas and n is undefined...
                        if (isfinite(m->tris[t].n) == false) continue;
                        //if (fabs(m->tris[t].n[0]) != fabs(m->tris[t].n[0])) continue;
                        //if (fabs(m->tris[t].n[1]) != fabs(m->tris[t].n[1])) continue;
                        //if (fabs(m->tris[t].n[2]) != fabs(m->tris[t].n[2])) continue;
                        //hash the normal
                        hash.add_point(mathtool::Point3d(m->tris[t].n));
                        hash.add_point(mathtool::Point3d(-m->tris[t].n));
                    }


                    for (uint v = 0; v < m->v_size; v++)
                    {
                        Vector3d vn;
                        for (auto fid : m->vertices[v].m_f)
                        {
                            vn = vn + m->tris[fid].n;
                        }
                        vn = vn.normalize();

                        if (isfinite(vn) == false) continue;

                        hash.add_point(mathtool::Point3d(vn));
                        hash.add_point(mathtool::Point3d(-vn));
                    }

                    all_dirs = hash.get_all_points();
                }

                //cluster directions if there are too many directions...
                if (all_dirs.size() > N) {
                    vector<uint> labels;
                    vector<mathtool::Point3d> centers;
                    normal_kmeans_clustering(N, all_dirs, labels, centers);
                    centers.swap(all_dirs);
                }
            }

            for (auto& dir : all_dirs)
            {
                N_dirs.push_back(Vector3d(dir.get()));
            }

            return N_dirs.size();
            //cout << "- Determine directions takes " << (getTime() - dir_start) / 1000 << " sec" << endl;
        }

        void normal_kmeans_clustering(int K, const vector<mathtool::Point3d> & all_dirs, vector<uint> & labels, vector<mathtool::Point3d> & centers)
        {
            int dir_size = all_dirs.size();

            cv::Mat samples(dir_size, 3, CV_32F);

            //init the samples
            for (int y = 0; y < dir_size; y++)
            {
                for (int z = 0; z < 3; z++)
                {
                    samples.at<float>(y, z) = all_dirs[y][z];
                }
            }

            //run k-means
            cv::Mat cvlabels;
            cv::Mat cvcenters;
            const int attempts = 10;
            double compactness =
                cv::kmeans(samples, K, cvlabels, cv::TermCriteria(CV_TERMCRIT_ITER | CV_TERMCRIT_EPS, 10000, 0.001), attempts, cv::KMEANS_PP_CENTERS, cvcenters);

            //parse output
            if (K > cvcenters.size().height) K = cvcenters.size().height;

            if (K == 0)
            {
                cout << "! [ClusterUnfolding] Error: normal_kmeans_clustering failed. No cluster found." << endl;
                return; //no clusters found?!
            }

            centers.reserve(K);
            for (int i = 0; i < K; i++)
            {
                float x = cvcenters.at<float>(i, 0);
                float y = cvcenters.at<float>(i, 1);
                float z = cvcenters.at<float>(i, 2);
                //cout << "center["<<i<<"]="<<Point3d(x, y, z) << endl;
                centers.push_back(mathtool::Point3d(x, y, z));
            }

            labels.reserve(dir_size);
            //vector< vector<uint> > clusters(K, vector<uint>() );
            for (int i = 0; i < dir_size; i++)
            {
                int id = cvlabels.at<int>(i, 0);
                //cout << i << "->" << id << endl;
                //clusters[id].push_back(i);
                labels.push_back(id);
            }

            //
            //for (int i = 0; i < K; i++) cout << "cluster " << i << " has " << clusters[i].size() << " elements" << endl;
            ///
        }

        //Using SPECTRAL CLUSTERING
        void ClusterUnfolding::spectral_score_clustering(int K, vector< vector<float> >& scores, vector<FaceCluster>& clusters)
        {
            cv::Mat idx;
            cv::Mat_<float> distMat;
            int tsize = scores.size();
            distMat.create(tsize, tsize);

            auto start = clock();
            for (int i = 0; i < tsize; i++)
            {
                for (int j = 0; j < tsize; j++)
                {
                    distMat(i, j) = 1 - scores[i][j];
                }
            }
            cout << "- [ClusterUnfolding] distMat computing time=" << (float)(clock() - start) / CLOCKS_PER_SEC << endl;


            //call open cv to cluster for us
            start = clock();
            cv::spectralClustering(distMat, idx, K, cv::LAPLACIAN_SHI_MALIK); //LAPLACIAN_SHI_MALIK, LAPLACIAN_NG_WEISS
            cout << "- [ClusterUnfolding] spectral clustering time=" << (float)(clock() - start) / CLOCKS_PER_SEC << endl;

            //convert the output to what we wanted.
            start = clock();
            vector<vector<uint> > loc_clusters(K, vector<uint>());
            int size = idx.size().height;

            for (int i = 0; i < size; i++) {
                uint label = idx.at<int>(i, 0);
                loc_clusters[label].push_back(i);
            }

            for (auto& cluster : loc_clusters)
            {
                FaceCluster fc;
                fc.fids = list<uint>(cluster.begin(), cluster.end());
                clusters.push_back(fc);
            }
            cout << "- [ClusterUnfolding] convert to desired format time=" << (float)(clock() - start) / CLOCKS_PER_SEC << endl;

            //start = clock();
            //save_score_to_image(loc_clusters);
            //cout << "save score to image time=" << (float)(clock() - start) / CLOCKS_PER_SEC << endl;
        }


        /// repair clusters so that each cluster contains a compact connected component
        void ClusterUnfolding::repair_clusters(model * m, vector<FaceCluster>& clusters, vector< vector<float> >& scores)
        {
            cout << "- [ClusterUnfolding] repair clusters" << endl;
            uint cluster_id = 0;

            //register the cluster id to each face
            for (auto& cluster : clusters)
            {
                for (auto fid : cluster.fids)
                {
                    m->tris[fid].cluster_id = cluster_id;
                }
                cluster_id++;
            }


            //get all connected components
            for (auto& cluster : clusters)
            {
                list<FaceCluster> ccs;
                get_connected_components(m, cluster, ccs);

                //cout << "cluster=" << m->tris[cluster.fids.front()].cluster_id << " ccs size=" << ccs.size()<< endl;
                //there is only one CC! great
                if (ccs.size() > 1) //there are multiple CCs
                {
                    //keep only the largest CC, which is the first in the list
                    list<uint> removed_fids;
                    for (auto ic = ++ccs.begin(); ic != ccs.end(); ic++) //for each cc
                    {
                        if(ic->fids.size()==1)
                        {
                            reassign_isolated_facet(m, ic->fids.front(), clusters, scores);
                        }
                        else //ic->fids.size()>1
                        {
                            for (auto fid : ic->fids) //for each face in cc
                            {
                                removed_fids.push_back(fid);
                            }
                        }
                    }

                    if (removed_fids.empty() == false) reassign_isolated_facets(m, removed_fids, clusters, scores);
                }
            }

            //find isolated and dangling facets
            for (uint fid = 0; fid < m->t_size; fid++)
            {
                int count = count_same_connected_cluster(m, fid);
                if (count == 0){//isolated
                    cout << "- [ClusterUnfolding] fid=" << fid << " cluster id="<< m->tris[fid].cluster_id<<" count == 0 should not happen" << endl;
                    assert(false);
                }

                if (count == 1)//dangling
                {
                    reassign_dangling_facet(m, fid, clusters, scores);
                }
            }//end for fid

            //recollect facet
            for (auto& cluster : clusters) cluster.fids.clear();
            for (uint fid = 0; fid < m->t_size; fid++)
            {
                clusters[m->tris[fid].cluster_id].fids.push_back(fid);
            }
            //done!!
        }

        /// get connected components of a give cluster
        /// results are stored in ccs
        void ClusterUnfolding::get_connected_components(model * m, FaceCluster& cluster, list<FaceCluster>& ccs) const
        {
            for (auto fid : cluster.fids)
            {
                m->tris[fid].parent_id = -1; //used later for collecting CC
            }

            double max_area = 0;
            FaceCluster * max_cc = NULL;
            for (auto fid : cluster.fids)
            {
                if (m->tris[fid].parent_id == -1)
                {
                    FaceCluster cc;
                    ccs.push_back(cc);
                    get_connected_component(m, fid, ccs.back());
                    double cc_area = cluster_area(m, ccs.back());
                    if (cc_area > max_area)
                    {
                        max_area = cc_area;
                        max_cc = &ccs.back();
                    }
                }
            }//end for fid

            //keep the largest component at the front
            if (max_cc != NULL && (&ccs.front()) != max_cc)
            {
                for (auto i = ccs.begin(); i != ccs.end(); i++)
                {
                    if (&(*i) == max_cc)
                    {
                        ccs.front().fids.swap(max_cc->fids);
                        break;
                    }
                }//end for i
            }//end fi
        }

        void ClusterUnfolding::get_connected_component(model * m, uint fid, FaceCluster& cc) const
        {
            uint my_cluster_id = m->tris[fid].cluster_id;
            m->tris[fid].parent_id = fid;

            list<uint> open;
            open.push_back(fid);
            while (open.empty() == false)
            {
                uint id = open.front();
                open.pop_front();
                cc.fids.push_back(id);

                auto neis = this->get_neighbor_triangles(m,id);
                for (auto nei : neis)
                {
                    if (m->tris[nei].parent_id == fid) continue;
                    if (m->tris[nei].cluster_id != my_cluster_id) continue;
                    m->tris[nei].parent_id = fid;
                    open.push_back(nei);
                }//end for nei
            }//end while
        }

        short ClusterUnfolding::count_same_connected_cluster(model * m, uint fid) const
        {
            short count = 0;
            auto neis = this->get_neighbor_triangles(m, fid);
            uint my_cluster_id = m->tris[fid].cluster_id;
            for (auto nei : neis)
            {
                if (m->tris[nei].cluster_id == my_cluster_id)
                    count++;
            }

            return count;
        }

        //
        // re-assign the facet to another cluster
        //
        void ClusterUnfolding::reassign_isolated_facet(model * m, uint fid, vector<FaceCluster>& clusters, const vector< vector<float> >& scores)
        {
            unordered_set<uint> adj_clusters;

            //collect adjacent cluster
            auto neis = this->get_neighbor_triangles(m, fid);
            for (auto nei : neis)
            {
                adj_clusters.insert(m->tris[nei].cluster_id);
            }

            //
            if (adj_clusters.size() == 1)
            {
                uint old = m->tris[fid].cluster_id;
                m->tris[fid].cluster_id = *adj_clusters.begin();
                cout << "- [ClusterUnfolding] reassign_isolated_facet fid=" << fid << " changed from " << old<<" to cluster = " << m->tris[fid].cluster_id << endl;
            }
            else //more than one possible facet
            {
                reassign_facet(m, fid, adj_clusters, clusters, scores);
            }
        }


        //
        // re-assign multiple connected isolated facets to another cluster
        // assume that all faces in fids are from the same cluster
        //
        void ClusterUnfolding::reassign_isolated_facets(model * m, list<uint>& fids, vector<FaceCluster>& clusters, const vector< vector<float> >& scores)
        {
            uint my_cid = m->tris[fids.front()].cluster_id;

            list<uint> open = fids;
            while (open.empty() == false) //loop through all facets until they are all handled
            {
                uint fid = open.front();
                open.pop_front();

                //collect adjacent cluster
                unordered_set<uint> adj_clusters;
                auto neis = this->get_neighbor_triangles(m, fid);
                for (auto nei : neis)
                {
                    if (m->tris[nei].cluster_id == my_cid) continue;
                    adj_clusters.insert(m->tris[nei].cluster_id);
                }

                if (adj_clusters.empty()) //not handled... wait until later
                {
                    open.push_back(fid);
                    continue;
                }

                //get a new cluster for this face
                int old_cid = m->tris[fid].cluster_id;
                reassign_facet(m, fid, adj_clusters, clusters, scores);
                if (old_cid == m->tris[fid].cluster_id)
                {
                    cout << "- [ClusterUnfolding] reassign_facet failed old_cid= " << old_cid << " my_cid= " << my_cid << endl;
                }
            }
        }

        void ClusterUnfolding::reassign_dangling_facet(model * m, uint fid, vector<FaceCluster>& clusters, const vector< vector<float> >& scores)
        {
            unordered_set<uint> adj_clusters;

            //collect adjacent cluster
            uint my_cid = m->tris[fid].cluster_id;

            auto neis = this->get_neighbor_triangles(m, fid);
            for (auto nei : neis)
            {
                adj_clusters.insert(m->tris[nei].cluster_id);
            }

            //
            if (adj_clusters.size() == 3)
            {
                reassign_facet(m, fid, adj_clusters, clusters, scores);
            }
            else //adj_clusters.size() == 2 (this is the only case)
            {
                assert(adj_clusters.size() == 2);
                adj_clusters.erase(my_cid);
                m->tris[fid].cluster_id = *adj_clusters.cbegin();
            }//end if

            //if this face is reassigned, check if the neighors become dangling
            if (m->tris[fid].cluster_id != my_cid)
            {
                //now check if the neighbors are also dangleing...
                for (auto nei : neis)
                {
                    if (m->tris[nei].cluster_id != my_cid) continue;
                    auto count = count_same_connected_cluster(m, nei);
                    if (count == 1)
                    {
                        reassign_dangling_facet(m, nei, clusters, scores); //handle this facets
                    }
                }//end for nei
            }

        }

        double ClusterUnfolding::avg_dist_to_cluster(model * m, uint fid, FaceCluster& cluster, const vector< vector<float> >& scores)
        {
            double total = 0;
            double total_area = 0;

            for (auto fid2 : cluster.fids)
            {
                double A = area(m, m->tris[fid2]);
                total += scores[fid][fid2] * A;
                total_area += A;
            }

            return total / total_area;
        }

        void ClusterUnfolding::reassign_facet
            (model * m, uint fid, unordered_set<uint> adj_clusters, vector<FaceCluster>& clusters, const vector< vector<float> >& scores)
        {
            uint best_cid = m->tris[fid].cluster_id;
            double min_error = DBL_MAX;
            for (uint cid : adj_clusters)
            {
                double error = avg_dist_to_cluster(m, fid, clusters[cid], scores);
                if (error < min_error)
                {
                    min_error = error;
                    best_cid = cid;
                }
            }//end for cid

            uint old = m->tris[fid].cluster_id;

            m->tris[fid].cluster_id = best_cid;

            cout << "- [ClusterUnfolding] reassign_facet fid=" << fid << " changed from " << old << " to cluster = " << m->tris[fid].cluster_id << endl;
        }

        double ClusterUnfolding::cluster_area(model * m, const FaceCluster& cluster) const
        {
            double total_area=0;
            for (auto fid : cluster.fids)
            {
                double A = area(m, m->tris[fid]);
                total_area += A;
            }
            return total_area;
        }

        //
        //
        // data serialization
        //
        //

        /// filename for serialization
        string ClusterUnfolding::get_serialization_filename() const
        {
          const auto config = this->m_unfolder->getConfig();
            auto m = this->m_unfolder->getModel();
            stringstream ss;
            ss.precision(2);
            ss << "." << this->m_unfolder->getModel()->name<<"_cluster_unfolding_score_";
            ss << m->v_size << "_" << m->t_size << "_" << m->e_size;
            ss << ".txt";

            if(this->m_unfolder->getConfig().binary_format)
      {
        return "." + m->name + ".bin";
      }

            if(config.scalar_score)
      {
        return "." + m->name + ".score";
      }

            return ss.str();
        }



        /// save score to file using filename from get_serialization_filename()
        void ClusterUnfolding::save_score_to_file() const
        {

            string filename = get_serialization_filename();
            ofstream fout;

            cout<<"saving score to "<<filename<<endl;

            const auto cfg = this->m_unfolder->getConfig();

            if(cfg.binary_format)
            {
              fout.open(filename.c_str(), std::ios::out | std::ios::binary);
            }else{
              fout.open(filename.c_str(), std::ios::out);
            }

            if (fout.good() == false)
            {
                cerr << "! [ClusterUnfolding] Error: save_score_to_file: Failed to open file " << filename << " to write." << endl;
                return;
            }

            fout.precision(10);

            //
            auto m = this->m_unfolder->getModel();

            if(cfg.scalar_score)
            {
              for(int i=0;i<m->t_size;++i)
                fout << this->m_score_vector[i] << endl;
            }
            else
            {
        // each weight is 4 bytes, should be smaller when number is large than 100 or number is float
        if(cfg.binary_format)
        {
          float runs = m_runs;
          float t_size = m->t_size;

          vector<int> ss(t_size);
          for(int i=0;i<t_size;++i)
            ss[i] = i;

//          std::random_shuffle(ss.begin(), ss.end());

          fout.write(reinterpret_cast<const char*>( &runs ), sizeof( float ));
          fout.write(reinterpret_cast<const char*>( &t_size ), sizeof( float ));
          for (int i = 0; i < m->t_size; i++)
          {
            int ii = ss[i];
            for (int j = i + 1; j < m->t_size; j++)
            {
              int jj = ss[j];
              fout.write(reinterpret_cast<const char*>( &this->m_score_matrix[ii][jj] ), sizeof( float ));
            }//end for j
          }//end for i
        }
        else
        {
          fout << this->m_runs << " ";
          fout << m->t_size << " ";
          for (int i = 0; i < m->t_size; i++)
          {
            for (int j = i + 1; j < m->t_size; j++)
            {
              fout << this->m_score_matrix[i][j] << " ";
            }//end for j
          }//end for i
        }
            }
        }

        ///load score from file using filename from get_serialization_filename()
        bool ClusterUnfolding::load_score_from_file()
        {
            string filename = get_serialization_filename();
            ifstream fin;

            const auto cfg = this->m_unfolder->getConfig();

            if(cfg.binary_format)
      {
        fin.open(filename.c_str(), std::ios::in | std::ios::binary);
      }else{
        fin.open(filename.c_str(), std::ios::in);
      }

            if (fin.good() == false)
            {
                cerr << "! - [ClusterUnfolding] Warning: load_score_from_file: Failed to open file " << filename << endl;
                return false;
            }

            auto m = this->m_unfolder->getModel();
            int t_size = 0;
            uint runs = 0;

            float buff;

            if(cfg.binary_format)
            {
              // in binary format all are float
              fin.read((char*)&buff, 4);
              runs = (int)buff;
              fin.read((char*)&buff, 4);
              t_size = (int)std::round(buff);
            }
            else
            {
              fin >> runs >> t_size;
            }
            if (m->t_size != t_size) {
              cout <<"t_size does not match! expected:"<<m->t_size<<" actual:"<<t_size<<" at "<<__FILE__<<" "<<__LINE__<<endl;
              return false; //for some reasons the size does not match...
            }

            this->m_runs += runs;
            cout<<"loaded runs = "<<runs<<endl;

            for (int i = 0; i < t_size; i++)
            {
                for (int j = i + 1; j < t_size; j++)
                {
					if(cfg.binary_format) {
						fin.read((char*)&buff, 4);
					} else {
						fin >> buff;
					}

                    m_score_matrix[j][i] = m_score_matrix[i][j] = buff;
                }//end for j
            }//end for i


            return true;
        }

        void ClusterUnfolding::dumpObj(const string& path, model * m) const
        {
            cout << "- [ClusterUnfolding] dumping obj to " << path << endl;
            ofstream out(path);

            vector<Vector3d> vs(m->t_size * 3);

            for (auto fid = 0; fid < m->v_size; ++fid) {
                const auto& p = m->vertices[fid].p;
                out << "v " << p[0] << " " << p[1] << " " << p[2] << endl;
            }

            for (auto fid = 0; fid < m->t_size; ++fid) {
                out << "f " << m->tris[fid].v[0] + 1 << " " << m->tris[fid].v[1] + 1 << " " << m->tris[fid].v[2] + 1 << "\n";
            }

            out.close();
        }
    }
}
