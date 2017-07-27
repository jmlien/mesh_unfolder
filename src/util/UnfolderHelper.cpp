/*
 * UnfolderHelper.cpp
 *
 *  Created on: Mar 4, 2016
 *      Author: zxi
 */

#include "UnfolderHelper.h"

#include <unordered_set>
#include <algorithm>
#include "unfolder.h"
#include "model.h"
#include "util/SVGWriter.h"

namespace masc {
    namespace unfolding {
        namespace util {

            UnfolderHelper::UnfolderHelper(Unfolder* unfolder) {
                this->m_unfolder = unfolder;
                this->m_num_children = this->computeNumChildren();

            }

            float UnfolderHelper::computeBorderCutsLength() {
                auto m = this->m_unfolder->getNet();
                //cout << "m edge size = " << m->e_size<< endl;
                SVGWriter writer(m,
                                 this->m_unfolder->getConfig());
                writer.Init();
                float len = 0.0f;
                const auto border_cut_eids = writer.GetBorderCuts();
                for(uint i : *border_cut_eids) {
                    len += m->edges[i].length;
                }
                return len;
            }

            float UnfolderHelper::computeVertexDistFromSameParent() {
                m_unfolder->unfoldTo(0.0);
                const auto m = this->m_unfolder->getNet();

                const auto & vertices = m->vertices;

                vertex v_default;

                unordered_map< uint, vector<uint> > v_share_parents;
                for(uint i = 0; i < vertices.size(); i++) {
                    uint p_id = vertices[i].parent_id;
                    //cout << "pid = " << p_id <<endl;
                    while(p_id != v_default.parent_id &&
                          vertices[p_id].parent_id != v_default.parent_id) {

                        p_id = vertices[p_id].parent_id;
                    }

                    if(p_id != v_default.parent_id) {
                        v_share_parents[p_id].push_back(i);
                    }

                }

                float dist = 0.0f;
                int n = 0;
                for(auto const& v : v_share_parents) {
                    const vector< uint >& v_ids = v.second;
                    if(v_ids.size() >1) {
                        for(int i = 0; i < v_ids.size()-1; i++){
                            for(int j = i+1; j < v_ids.size(); j++) {
                                dist += (vertices[v_ids[i]].p -
                                         vertices[v_ids[j]].p).norm();
                                n++;

                            }
                        }
                    }
                }
                return dist;
            }

            float UnfolderHelper::computeVertexDist() {

                const auto & org = this->m_unfolder->getOrg();
                m_unfolder->unfoldTo(1.0);
                const auto & unfolded = this->m_unfolder->getUnfolded();
                cout << unfolded.size() << endl;
                vector<Vector3d> v_org(org.size()*3);
                vector<Vector3d> v_unfolded(unfolded.size()*3);

                for (int i = 0; i < org.size(); i++) {
                    for (int j = 0; j < 3; j++) {

                        v_org.push_back(org[i][j].second);
                        v_unfolded.push_back(unfolded[i][j].second);
                    }
                }

                float dist = 0;
                int n = 0;
                for (int i = 0; i < v_org.size(); i++) {
                    for ( int j = i+1; j < v_org.size(); j++) {
                        float org_dist = (v_org[i]-v_org[j]).norm();
                        float unfolded_dist = (v_unfolded[i]-v_unfolded[j]).norm();
                        //cout << org_dist << " " << unfolded_dist  << endl;
                        if (unfolded_dist < 0.01 && org_dist - unfolded_dist> 0.01) {
                            dist += unfolded_dist;
                            n++;
                        }
                    }
                }
                return dist/(float)n;
            }

// compute pairwise face geodesic distance in the unfolding
            GEO_DIST_GRAPH UnfolderHelper::computeGeoDist() {
                GEO_DIST_GRAPH dist;
                const auto m = this->m_unfolder->getModel();
                for (int i = 0; i < m->t_size; ++i)
                    for (int j = i + 1; j < m->t_size; ++j) {
                        dist[i][j] = this->computeGeoDist(i, j);
                    }
                return dist;
            }

            int UnfolderHelper::computeDiameter() {

                int diameter = 0;
                for (int i = 0; i < m_num_children.size(); i++) {
                    if (m_num_children[i] > 0) continue;
                    for (int j = i+1; j < m_num_children.size(); j++ ) {
                        if (m_num_children[j] == 0) continue;
                        int dist = computeDist(i,j);
                        if (dist > diameter) diameter = dist;
                    }
                }

                return diameter;
            }


// compute pairwise face geodesic distance for a given connected component
            std::vector<float> UnfolderHelper::computeGeoDist(const int root,
                                                              const std::vector<unsigned int>& cc) {

//  std::cout << "UnfolderHelper::computeGeoDist root = " << root << std::endl;

                std::vector<float> dist(cc.size(), 0.0f);
                const auto m = this->m_unfolder->getModel();
                for (int i = 0; i < cc.size(); ++i) {
                    int fid1 = std::min((int) cc[i], root);
                    int fid2 = std::max((int) cc[i], root);
                    if (fid1 == fid2)
                        continue;
                    dist[i] = this->computeGeoDist(fid1, fid2);
                }
                return dist;
            }

// compute dist of two faces
            int UnfolderHelper::computeDist(const int fid1, const int fid2) const {
                if (fid1 == fid2)
                    return 0;

                //visited parents
                std::unordered_set<int> p;
                int pid1 = fid1;
                int pid2 = fid2;

                // dist from fid1 to a parent <pid, dist>
                unordered_map<int, int> dist1;
                // dist from fid2 to a parent <pid, dist>
                unordered_map<int, int> dist2;

                p.insert(pid1);
                p.insert(pid2);

                dist1[pid1] = 0;
                dist2[pid2] = 0;

                const auto m = this->m_unfolder->getModel();

                int common_ancestor  = -1;

                while (true) {
                    int new_pid1 = m_unfolder->getParentFaceId(pid1);
                    if (new_pid1 >= 0) {
                        dist1[new_pid1] = dist1[pid1] + 1;
                        if (p.count(new_pid1)) {
                            common_ancestor = new_pid1;
                            break;
                        }
                        pid1 = new_pid1;
                        p.insert(pid1);
                    }
                    int new_pid2 = m_unfolder->getParentFaceId(pid2);
                    if (new_pid2 >= 0) {
                        dist2[new_pid2] = dist2[pid2] + 1;

                        if (p.count(new_pid2)) {
                            common_ancestor = new_pid2;
                            break;
                        }

                        pid2 = new_pid2;
                        p.insert(pid2);

                    }
                }
                return dist1[common_ancestor] + dist2[common_ancestor];
            }
            float UnfolderHelper::computeGeoDist(const int fid1, const int fid2) const{
//  std::cout << "UnfolderHelper::computeGeoDist " << fid1 << " " << fid2
//      << std::endl;

                if (fid1 == fid2)
                    return 0.0f;

                // vistited parents
                std::unordered_set<int> p;
                int pid1 = fid1;
                int pid2 = fid2;
                // dist from fid1 to a parent <pid, dist>
                unordered_map<int, float> dist1;
                // dist from fid2 to a parent <pid, dist>
                unordered_map<int, float> dist2;

                p.insert(pid1);
                p.insert(pid2);

                dist1[pid1] = 0.0f;
                dist2[pid2] = 0.0f;

                const auto m = this->m_unfolder->getModel();

                int common_ancestor = -1;

                while (true) {

                    int new_pid1 = m_unfolder->getParentFaceId(pid1);

                    if (new_pid1 >= 0) {
                        dist1[new_pid1] = dist1[pid1]
                            + (m->tris[new_pid1].center - m->tris[pid1].center).norm();

                        if (p.count(new_pid1)) {
                            common_ancestor = new_pid1;
                            break;
                        }

                        pid1 = new_pid1;
                        p.insert(pid1);
                    }

                    int new_pid2 = m_unfolder->getParentFaceId(pid2);
                    if (new_pid2 >= 0) {
                        dist2[new_pid2] = dist2[pid2]
                            + (m->tris[new_pid2].center - m->tris[pid2].center).norm();

                        if (p.count(new_pid2)) {
                            common_ancestor = new_pid2;
                            break;
                        }

                        pid2 = new_pid2;
                        p.insert(pid2);
                    }

                }

                float dist = dist1[common_ancestor] + dist2[common_ancestor];

//  std::cout << "dist = " << dist << std::endl;

                return dist;
            }

            int UnfolderHelper::computeNumLeafNodes() {


                int count = 0;
                for (auto& c : m_num_children) {
                    if (c == 0) count++;
                    //cout << c << std::endl;
                }

                return count;
            }


            std::vector<int> UnfolderHelper::computeBranchLength() {

                std::vector<int> b_len;
                b_len.resize(m_num_children.size());
                for(int i = 0; i < m_num_children.size(); i++) {
                    if(m_num_children[i] > 0) {
                        //not a leaf
                        b_len[i] = -1;
                        continue;
                    }

                    int new_pid = i;
                    while (true) {
                        new_pid = m_unfolder->getParentFaceId(new_pid);

                        if (new_pid < 0) break;

                        b_len[i] += 1;

                        if(m_num_children[new_pid] > 1) break;
                    }
                }
                return b_len;
            }

            std::vector<float> UnfolderHelper::computeBranchGeoLength() {

                const auto & m = this->m_unfolder->getModel();

                std::vector<float> b_len;
                b_len.resize(m->t_size);
                for(int i = 0; i < m->t_size; i++) {
                    if(m_num_children[i] > 0) {
                        //not a leaf
                        b_len[i] = -1.0f;
                        continue;
                    }
                    int pid = i;

                    while (true) {
                        int new_pid = m_unfolder->getParentFaceId(pid);

                        if (new_pid < 0) break;

                        b_len[i] += (m->tris[new_pid].center - m->tris[pid].center).norm();

                        if (m_num_children[new_pid] > 1) break;

                        pid = new_pid;
                    }
                }
                return b_len;
            }

            std::vector<int> UnfolderHelper::computeNumChildren() {

                std::vector<int> c;

                const auto m = this->m_unfolder->getModel();

                //cout << "t_size =!========= " << m->t_size <<endl;
                c.resize(m->t_size);
                for(int i = 0; i < m->t_size; ++i) {
                    int new_pid = m_unfolder->getParentFaceId(i);
                    //cout <<"new_pid = " << new_pid <<endl;
                    if (new_pid >= 0) {
                        c[new_pid] += 1;
                    }
                }

                return c;
            }

            float UnfolderHelper::getTotalCutLength() {

                float m_total_edge_length = 0.0f;

                const auto m = this->m_unfolder->getModel();
                // add cache for each model
                for (auto i = 0; i < m->e_size; i++) {
                    const auto& edge = m->edges[i];
                    const auto fid1 = edge.fid[0];
                    const auto fid2 = edge.fid[1];

                    const auto& v1 = m->vertices[edge.vid[0]];
                    const auto& v2 = m->vertices[edge.vid[1]];

                    const auto edge_len = (float) ((v2.p - v1.p).norm());

                    m_total_edge_length += edge_len;
                }

                float total_selected_edge_length = 0.0f;

                for (const auto& eid : m_unfolder->getFoldEdges()) {
                    const auto& e = m->edges[eid];
                    const auto& v1 = m->vertices[e.vid[0]];
                    const auto& v2 = m->vertices[e.vid[1]];

                    const auto edge_len = (float) ((v2.p - v1.p).norm());

                    total_selected_edge_length += edge_len;
                }

                return m_total_edge_length - total_selected_edge_length;
            }

            float UnfolderHelper::getTotalCutLengthNormalized() {

                // float m_total_edge_length = 0.0f;

                const auto m = this->m_unfolder->getModel();
                vector<float> edge_lens;
                edge_lens.reserve(m->e_size);
                // add cache for each model
                for (auto i = 0; i < m->e_size; i++) {
                    edge_lens.push_back(m->edges[i].length);
                    //m_total_edge_length += m->edges[i].length;
                }

                std::sort (edge_lens.begin(), edge_lens.end());
                float total_selected_edge_length = 0.0f;

                int edge_count = 0;
                for (const auto& eid : m_unfolder->getFoldEdges()) {
                    edge_count++;
                    total_selected_edge_length += m->edges[eid].length;
                }
                float min_cut_length = 0.0;
                float max_cut_length = 0.0;
                for (int i = 0; i < edge_count; i++) {
                    min_cut_length += edge_lens[i];
                    max_cut_length += edge_lens[edge_lens.size()-i-1];
                }
                //cout << "cut length = " << total_selected_edge_length << endl;

                //cout << "range = " << max_cut_length - min_cut_length << endl;
                return (total_selected_edge_length - min_cut_length)
                    / (max_cut_length - min_cut_length);
            }
            
           void UnfolderHelper::dumpStats(std::ostream & out){
                
                this->m_unfolder->rebuildModel();
                
                int num_leaf_nodes=0;
                int num_1_degree_nodes=0;
                int num_2_degree_nodes=0;
                int num_3_degree_nodes=0;
                int total_degrees=0;
                
                for (int i = 0; i < this->m_num_children.size(); i++) {
                    
                    switch(this->m_num_children[i]) {
                        case 0: num_leaf_nodes++; break;
                        case 1: num_1_degree_nodes++; break;
                        case 2: num_2_degree_nodes++; break;
                        case 3: num_3_degree_nodes++; break;
                    }
                    total_degrees += this->m_num_children[i];
                    
                }
                out << "w_cutlength " << this->getTotalCutLengthNormalized() << endl;
                out << "w_num_leaf_nodes    " <<  (double) num_leaf_nodes / (double)total_degrees << endl;
                out << "w_num_1_degree_nodes    " << (double) num_1_degree_nodes / (double)total_degrees << endl;
                out << "w_num_2_degree_nodes    " << (double) num_2_degree_nodes / (double)total_degrees << endl;
                out << "w_border_cuts   " << (double)this->computeBorderCutsLength() / (double)this->getTotalCutLength() << endl;
                out << "w_hullarea   " << this->m_unfolder->getModel()->surface_area / this->m_unfolder->getHullArea() << endl;
                
            }
            
            
            void UnfolderHelper::dumpStats(const std::string& path) {
                std::cout << "- dumping stats to "<< path << std::endl;

                std::ofstream out(path);
                this->dumpStats(out);
            }
        } /* namespace util */
    } /* namespace unfolding */
} /* namespace masc */
