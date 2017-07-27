
/*
* LPUnfolding.cpp
*
*  Created on: Feb 29, 2015
*      Author: jmlien
*/

#include "LPUnfolding.h"


#include <glpk.h> //linear programming solver
#include "libga/util/StreamHelper.h"

#include <deque>

namespace masc 
{
	namespace unfolding 
	{

		const bool debug = false;
		const bool dump_cluster_obj = true;
		
		//LP functions
		void callback(glp_tree *tree, void *info);
		int add_lazy_constraints(LPUnfolding* unfolder, glp_prob * lp);


		LPUnfolding::LPUnfolding(Unfolder* unfolder)
		{
			Config mycfg = unfolder->getConfig();
			mycfg.record_overlap = true;

			this->m_unfolder = new Unfolder(unfolder->getModel(), mycfg);
			this->m_unfolder_from_caller = unfolder;
			assert(this->m_unfolder);
			m_visited = new unsigned char[unfolder->getModel()->t_size];
			assert(this->m_visited);
		}

		LPUnfolding::~LPUnfolding()
		{
			delete this->m_unfolder;
			delete this->m_visited;
		}


		void LPUnfolding::init()
		{
			srand48(time(NULL));
		}

		bool LPUnfolding::setup(const string& filename)
		{
			ifstream in;
			in.open(filename, ios_base::in);
			if (!in.good())
			{
				cerr << "!Error in LPUnfolding::setup. Can not open file = " << filename << endl;
				return false;
			}

			// invoke implementation
			auto flag = this->setup(in);

			in.close();

			return flag;
		}

		bool LPUnfolding::setup(istream& in)
		{
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
				}
				else
				{
					break;
				}
			}

			this->init();

			//srand48(time(nullptr));

			this->print(cout);


			return true;
		}

		bool LPUnfolding::setup(vector<string>& tokens)
		{
			//this->setSeed(time(nullptr));

			for (int i = 0; i<tokens.size(); ++i)
			{
				const auto& token = tokens[i];
			}
			return true;
		}

		//////////////////////////////////////////
		// run n unfoldings and cluster
		//////////////////////////////////////////

		void LPUnfolding::run()
		{
			auto start_time = clock();
			vector<LP_constraints> constaints;
			

			model * mesh = this->m_unfolder->getModel();
			this->GenerateConstraints(mesh, constaints);

			while (true)
			{
				vector<float> solution;
				bool solved = SolveLP(mesh, constaints, solution);

				if (solved)
				{
					cout << "- Optimal solution found. ";
					int added_constaints = this->AddLazyConstraints(mesh, solution, constaints);
					/*
					cout << "- Optimal solution found: ";
					int id = 0;
					double sum = 0;

					cout << "Weight = ";
					for (auto x : solution) cout << x << " ";
					cout << endl;

					for (auto x : solution)
					{
					sum += x;
					cout << "[" << id++ << "]="<<x << endl;
					}
					cout << endl;
					cout << "# of connecting edge = " << sum << "/" << mesh->e_size << endl;;
					cout << "total tri size = " << mesh->t_size << endl;

					//triangle constaints (not all edges can be cut edge)
					for (int fid = 0; fid < mesh->t_size; fid++)
					{
					auto & tri = mesh->tris[fid];
					int cut_count = 0;
					for (short d = 0; d < 3; d++)
					{
					cut_count+=solution[tri.e[d]];
					}
					cout << "f[" << fid << "] has " << 3 - cut_count << " cuts" << endl;
					}

					//vertex constaints (at least one cut for regular vertex
					//and two cuts for hyperbolic vertex)
					for (int vid = 0; vid < mesh->v_size; vid++)
					{
					auto & vertex = mesh->vertices[vid];
					int cut_count = 0;
					for (auto eid : vertex.m_e)
					{
					cut_count += solution[eid];
					}
					cout << "v[" << vid << "] is "<<(vertex.hyperbolic?"hypo":"regular")<<" has " << vertex.m_e.size()- cut_count << " cuts" << endl;

					}
					*/

					if (added_constaints == 0)
					{
						//tell the unfolder from the caller to unfold...
						this->m_unfolder_from_caller->buildFromWeights(solution);
						this->m_unfolder_from_caller->rebuildModel();
						cout << "No Overlaps!!" << endl;
						break;
					}
					else
					{
						cout << "\n====================================================" << endl;
						cout << "Find " << added_constaints << " overlaping pairs. Keep working..." << endl;
						cout << "constaints size=" << constaints.size() << endl;
						cout << "====================================================" << endl;
					}
				}
				else
				{
					cout << "- Optimal solution NOT found" << endl;
					break;
				}
			}
		}

		//////////////////////////////////////////

		void LPUnfolding::print(ostream& out) const
		{

		}


		//////////////////////////////////////////
		// LP
		//////////////////////////////////////////
		
#if 0
		void callback_test(glp_tree *tree, void *info)
		{
			cout << "-----------------------------------------------" << endl;

			int a_cnt=0,n_cnt=0,t_cnt=0;
			glp_ios_tree_size(tree, &a_cnt, &n_cnt, &t_cnt);			cout << "tree size, active nodes= " << a_cnt << " all=" << n_cnt << " total=" << t_cnt << endl;

			Unfolder* unfolder = dynamic_cast<Unfolder*>((Unfolder*)info);
			assert(unfolder);
			glp_prob * prob = glp_ios_get_prob(tree);
			int curr_node = glp_ios_curr_node(tree);
			int best_node = glp_ios_best_node(tree);

			cout << "unfolder edge size=" << unfolder->getModel()->e_size << endl;
			cout << "curr_node=" << curr_node << " best_node=" << best_node << " prob=" << prob << endl;

			//
			if (curr_node != 0)
			{

				cout << "cur level = " << glp_ios_node_level(tree, curr_node) << endl;
				cout << "cur bound = " << glp_ios_node_bound(tree, curr_node) << endl;
				cout << "row size=" << glp_get_num_rows(prob) << endl;;
				////cout << "node data = " << glp_ios_node_data(tree, curr_node) << endl;;
				cout << "node data = " << glp_ios_node_data(tree, curr_node) << endl;
				//for (int i = 1; i <= glp_get_num_rows(prob); i++)
				//{
				//	glp_attr attr;
				//	glp_ios_row_attr(tree, i, &attr);
				//	cout << "row[" << i << "] level=" << attr.level << " origin=" << attr.origin << " klass=" << attr.klass << endl;
				//}

				//
				/* access subproblem application-specific data */
			}

			auto reason = glp_ios_reason(tree);

			if (reason == GLP_IPREPRO) cout << "Request for preprocessing" << endl;
			else if (reason == GLP_IROWGEN){
				int glp_prim_stat = glp_mip_status(prob);

				if (glp_prim_stat == GLP_FEAS)
				{
					cout << "add one or more \"lazy\" constraints" << endl;

					//int nrow = glp_add_rows(prob, 1);
					//int total_var = glp_get_num_cols(prob);
					//
					//int id = (total_var + 1)*drand48();
					//if (id < 1) id = 1;
					//if (id>total_var) id = total_var;

					//glp_set_row_bnds(prob, nrow, GLP_FX, 0, 0);
					//int index[2] = { 0, id };
					//double value[2] = { 0, 1 };
					//glp_set_mat_row(prob, nrow, 1, index, value);
				}
			}
			else if (reason == GLP_IHEUR) cout << "Request for heuristic solution" << endl; //glp_ios_heur_sol(tree, an_array); break;
			else if (reason == GLP_ICUTGEN)
			{
				cout << "Request for cut generation" << endl;
			}
			else if (reason == GLP_IBRANCH){
				cout << "Request for branch" << endl;
				//glp_ios_branch_upon(tree, 1, GLP_NO_BRNCH);
			}
			else if (reason == GLP_ISELECT) cout << "Request for subproblem selection" << endl;
			else if (reason == GLP_IBINGO) cout << "Bingo" << endl;



			int glp_prim_stat = glp_mip_status(prob);

			switch (glp_prim_stat)
			{
			case GLP_OPT: cout << "solution is optimal;" << endl; break;
			case GLP_FEAS: cout << "solution is feasible;" << endl; break;
			case GLP_INFEAS: cout << "solution is infeasible;" << endl; break;
			case GLP_NOFEAS: cout << "problem has no feasible solution;" << endl; break;
			case GLP_UNBND: cout << "problem has unbounded solution;" << endl; break;
			case GLP_UNDEF: cout << "solution is undefined." << endl; break;
			}


			//vector<double> current_solution;
			//cout << "current solution:";
			//for (int i = 1; i <= unfolder->getModel()->e_size; i++)
			//{
			//	double x = glp_mip_col_val(prob, i);
			//	current_solution.push_back(x);
			//	cout << x << ",";
			//}
			//cout << endl;
			
			//select subproblem to continue the search
			//glp_ios_select_node(tree,0);


		}
#endif

		bool LPUnfolding::SolveLP(model * m, vector<LP_constraints>& constaints, vector<float>& solution)
		{
			int constraint_size = constaints.size();
			int variable_size = m->e_size;

			glp_prob * lp = glp_create_prob();
			assert(lp);
			glp_set_prob_name(lp, "lp");
			glp_set_obj_dir(lp, GLP_MIN);
			glp_add_rows(lp, constraint_size);
			glp_add_cols(lp, variable_size);

			//init rows
			int iaja_size = 0;
			int row_id = 1;

			char tmp[64];
			for (auto & c : constaints)
			{
				sprintf(tmp, "r%08d", row_id);
				glp_set_row_name(lp, row_id, tmp);
				glp_set_row_bnds(lp, row_id, c.type, c.lower_bound, c.upper_bound);
				iaja_size += c.eids.size();
				row_id++;
			}

			//init cols
			for (int i = 1; i <= variable_size; i++)
			{
				char tmp[64];
				sprintf(tmp, "s%08d", i);
				glp_set_col_name(lp, i, tmp);
				//MIP
				glp_set_col_kind(lp, i, GLP_BV); 
				//MIP
				{
					auto & e = m->edges[i - 1];
					const auto  & vec = m->vertices[e.vid[0]].p - m->vertices[e.vid[1]].p;
					glp_set_obj_coef(lp, i, vec.norm());
				}

			}//end i

			//init ia, ja, and ar
			int * ia = new int[1 + iaja_size];
			int * ja = new int[1 + iaja_size];
			double * ar = new double[1 + iaja_size];
			assert(ia && ja && ar);

			int ia_id = 1;
			row_id = 1;
			for (auto & c : constaints)
			{
				for (auto eid : c.eids)
				{
					ia[ia_id] = row_id;
					ja[ia_id] = eid + 1;
					ar[ia_id] = 1;
					ia_id++;
				}//end for j

				row_id++;
			}//end for i


			glp_load_matrix(lp, iaja_size, ia, ja, ar);

			//assert(glp_simplex(lp, NULL) == 0);
			//assert(glp_get_status(lp) == GLP_OPT);

			glp_iocp parm;
			glp_init_iocp(&parm);
			parm.pp_tech = GLP_PP_ALL;
			parm.presolve = GLP_ON;
			parm.clq_cuts = GLP_ON;
			parm.binarize = GLP_ON;
			//parm.cb_func = callback;
			parm.cb_info = this;

			parm.tm_lim = 600000; //600 sec
			int err = glp_intopt(lp, &parm);
			cout << "err=" << err << endl;

			double z = glp_mip_obj_val(lp);
			cout << "objective value=" << z << endl;


			//get mip status
			int glp_prim_stat = glp_mip_status(lp);

			switch (glp_prim_stat)
			{
			case GLP_OPT: cout << "solution is optimal;" << endl; break;
			case GLP_FEAS: cout << "solution is feasible;" << endl; break;
			case GLP_INFEAS: cout << "solution is infeasible;" << endl; break;
			case GLP_NOFEAS: cout << "problem has no feasible solution;" << endl; break;
			case GLP_UNBND: cout << "problem has unbounded solution;" << endl; break;
			case GLP_UNDEF: cout << "solution is undefined." << endl; break;
			}

			bool solution_found = glp_prim_stat == GLP_OPT || glp_prim_stat == GLP_FEAS;

			if (solution_found)
			{
				for (int i = 1; i <= variable_size; i++)
				{
					double x = glp_mip_col_val(lp, i);
					solution.push_back(x);
				}
			}

			glp_delete_prob(lp);
			delete[] ia;
			delete[] ja;
			delete[] ar;

			return solution_found;
		}

		void LPUnfolding::GenerateConstraints(model * m, vector<LP_constraints>& constaints)
		{
			//  The constraints would be, for each triangle, the sum of edge cut is $ <= 2$,
			//  for each vertex, the sum is $ >= 2$ for hyperbolic vertex or $ >= 1$ otherwise.
			//	We also constrain that the number of the cut edges is exactly $F - 1$ to ensure that the net remains in a single connected component.
			//	The objective function is used to minimize the total cut length.


			//triangle constaints (not all edges can be cut edge)
			for (int fid = 0; fid < m->t_size; fid++)
			{
				auto & tri = m->tris[fid];
				LP_constraints c;
				for (short d = 0; d < 3; d++)
				{
					c.eids.push_back(tri.e[d]);
				}
				c.type = GLP_DB;
				c.lower_bound = 0;
				c.upper_bound = 2;
				constaints.push_back(c);
			}

			//vertex constaints (at least one cut for regular vertex 
			//and two cuts for hyperbolic vertex)
			for (int vid = 0; vid < m->v_size; vid++)
			{
				auto & vertex = m->vertices[vid];
				if (vertex.m_e.empty()) continue;

				LP_constraints c;
				c.eids= vector<uint>(vertex.m_e.begin(), vertex.m_e.end());
				c.type = GLP_LO;
				c.lower_bound = (vertex.hyperbolic)?2:1; //?
				constaints.push_back(c);

				//cout << "v[" << vid << "] ";
				//for (auto eid : c.eids) cout << eid << ", ";
				//cout << " >= " << c.lower_bound;
				//cout << endl;
			}

			//edge constaints
			LP_constraints c;
			for (int eid = 0; eid < m->e_size; eid++)
			{
				c.eids.push_back(eid);
			}
			c.type = GLP_FX;
			c.lower_bound = c.upper_bound = m->e_size - m->t_size + 1; //?

			constaints.push_back(c);
		}


		int LPUnfolding::AddLazyConstraints(model * m, vector<float>& weights, vector<LP_constraints>& constaints)
		{
			//----------------------------------------------------
			// need to make sure that one cc of faces
			//
			list< list<uint> > ccs;
			int ccsize = find_face_ccs(weights, ccs);

			if (ccsize > 0)
			{
				for (const auto & cc : ccs)
				{
					LP_constraints c;
					c.eids = vector<uint>(cc.begin(), cc.end());
					c.upper_bound = cc.size()-1;
					c.type = GLP_UP;
					constaints.push_back(c);
				}
				return ccsize; //well this weight cannot produce a mst so we are done here
			}

			return 0;

			//----------------------------------------------------
			//OK, now we have one cc....add constraints to avoid overlapping
			list< list<uint> > loops;
			find_conflict_loops(weights, loops);
			auto loop_size = loops.size();

			for (const auto & loop : loops)
			{
				LP_constraints c;
				c.eids=vector<uint>(loop.begin(),loop.end());
				c.lower_bound = 1;
				c.type = GLP_LO;
				constaints.push_back(c);
			}//end for loop

			return loop_size;
		}


		//this function add additional overlapping constraints to the given program.
		int add_lazy_constraints(LPUnfolding* unfolder, glp_prob * lp)
		{
			int total_var = glp_get_num_cols(lp);
			vector<float> weights;
			for (int i = 1; i <= total_var; i++)
			{
				float x = (float)glp_mip_col_val(lp, i);
				//float x = (float)glp_get_col_prim(lp, i);
				if (x < 0.1) x = 0;
				else x = 1;
				weights.push_back(x);
			}


			// need to make sure that one cc of faces
			//
			list< list<uint> > ccs;
			int ccsize = unfolder->find_face_ccs(weights, ccs);
			int nrow = glp_add_rows(lp, ccsize);

			cout << "nrow=" << nrow << " ";

			if (ccsize > 0)
			{
				for (const auto & cc : ccs)
				{
					int csize = cc.size();
					cout << csize << " ";
					glp_set_row_bnds(lp, nrow, GLP_UP, 0, csize-1);

					int * index = new int[csize + 1];
					double * value = new double[csize + 1];
					int id = 1;
					for (const auto & eid : cc)
					{
						index[id] = eid + 1;
						value[id] = 1;
						id++;
					}
					glp_set_mat_row(lp, nrow, csize, index, value);
					delete[] index;
					delete[] value;
					nrow++;
				}
				cout << endl;
				return ccsize; //well this weight cannot produce a mst so we are done here
			}
			return 0;

			//find loops
			list< list<uint> > loops;
			unfolder->find_conflict_loops(weights, loops);

			//no conflicts!!
			if (loops.empty()) return 0;

			//for each loop, we create a constraint, such that the look can be broken
			//that is, one of the edges in the loop must be cut

			//try add constraint directly (or we can add the constraint to the cut pool)
			auto loop_size = loops.size();
			nrow = glp_add_rows(lp, loop_size);
			for (const auto & loop : loops)
			{
				glp_set_row_bnds(lp, nrow, GLP_LO, 1, 0);
				int csize = loop.size();
				int * index = new int[csize + 1];
				double * value = new double[csize + 1];
				int id = 1;
				for (const auto & eid : loop)
				{
					index[id] = eid + 1;
					value[id] = 1;
					id++;
				}
				glp_set_mat_row(lp, nrow, csize, index, value);

				nrow++;
			}//end for loop

			return loop_size;
		}

		//generate lazy constraints here
		void callback(glp_tree *tree, void *info)
		{
			//cout << "-----------------------------------------------" << endl;

			LPUnfolding* unfolder = dynamic_cast<LPUnfolding*>((LPUnfolding*)info);
			assert(unfolder);

			glp_prob * lp = glp_ios_get_prob(tree);
			int curr_node = glp_ios_curr_node(tree);
			int best_node = glp_ios_best_node(tree);

			auto reason = glp_ios_reason(tree);

			//if (reason == GLP_IPREPRO) cout << "Request for preprocessing" << endl;
			//else if (reason == GLP_IROWGEN){ cout << "add one or more \"lazy\" constraints" << endl;}
			//else if (reason == GLP_IHEUR) cout << "Request for heuristic solution" << endl; //glp_ios_heur_sol(tree, an_array); break;
			//else if (reason == GLP_ICUTGEN){ cout << "Request for cut generation" << endl; }
			//else if (reason == GLP_IBRANCH){cout << "Request for branch" << endl;}
			//else if (reason == GLP_ISELECT) cout << "Request for subproblem selection" << endl;
			//else if (reason == GLP_IBINGO) cout << "Bingo" << endl;

			if (reason == GLP_IPREPRO)
			{
				//solve the problem to get a feasible solution and add additional constraints
				
			}
			else if (reason == GLP_IROWGEN)
			{
				int glp_prim_stat = glp_mip_status(lp);

				if (glp_prim_stat == GLP_FEAS) //make sure that the solution is feable
				{
					int a_cnt = 0, n_cnt = 0, t_cnt = 0;
					glp_ios_tree_size(tree, &a_cnt, &n_cnt, &t_cnt);					cout << "tree size, active nodes= " << a_cnt << " all=" << n_cnt << " total=" << t_cnt << " curr_node=" << curr_node << endl;

					//cout << "add one or more \"lazy\" constraints" << endl;
					int size=add_lazy_constraints(unfolder, lp);
					cout << "add " << size<<" \"lazy\" constraints" << endl;
				}//end if (glp_prim_stat == GLP_FEAS)
			}//end if reason == GLP_IROWGEN
		}//end callback(glp_tree *tree, void *info)

		//////////////////////////////////////////
		// find connected components in the mesh
		//////////////////////////////////////////

		//note that each cc is composed of a list of edge ids; not facet ids
		int LPUnfolding::find_face_ccs(const vector<float>& weight, list< list<uint> >& ccs)
		{
			//in weight: 1 mean cut and 0 mean no cut
			auto mesh = this->m_unfolder->getModel();
			memset(m_visited, 0, sizeof(unsigned char)* mesh->t_size);

			list< list<uint> > face_ccs;
			int ccsize = 0;
			for (uint fid = 0; fid < mesh->t_size; fid++)
			{
				if (m_visited[fid]) continue; //visited
				face_ccs.push_back(list<uint>());
				find_face_cc(fid, weight, face_ccs.back());
				ccsize++;
			}

			//assign cc id to each face
			uint ccid = 0;
			for (auto& cc : face_ccs)
			{
				for (uint fid : cc)
				{
					mesh->tris[fid].cluster_id = ccid;
				}
				ccid++;
			}

			//identify boundary edges of each cc
			for (auto& cc : face_ccs)
			{
				ccs.push_back(list<uint>());
				list<uint> & eids = ccs.back();
				for (uint fid : cc)
				{
					auto & tri = mesh->tris[fid];

					for (short d = 0; d < 3; d++)
					{
						uint eid = tri.e[d];
						if (weight[eid] <SMALLNUMBER) continue; //not cut edge
						uint ofid = get_neighbor_triangle(mesh, fid, d);
						if (mesh->tris[ofid].cluster_id == tri.cluster_id) continue; //same cluster...
						eids.push_back(eid);
					}//end for d
				}
			}

			return ccsize;
		}

		//get a single face cc
		void LPUnfolding::find_face_cc(uint start, const vector<float>& weight, list<uint> & cc)
		{
			//in weight: 1 means cut and 0 means no cut
			auto mesh = this->m_unfolder->getModel();
			std::deque<uint> open;
			m_visited[start] = true;
			open.push_back(start);

			while (open.empty() == false)
			{
				uint fid = open.front();
				open.pop_front();
				cc.push_back(fid);
				auto & tri = mesh->tris[fid];

				for (short d = 0; d < 3; d++)
				{
					uint eid = tri.e[d];
					if (weight[eid] >=0.999) continue; //cut edge
					uint ofid = get_neighbor_triangle(mesh, fid, d);
					if (m_visited[ofid]) continue; //visited
					m_visited[ofid] = true;
					open.push_back(ofid);
				}//end for d
			}//end while
		}

		//////////////////////////////////////////
		// find overlapping loops in the unfolding
		//////////////////////////////////////////

		//note that each loop is composed of a list of edge ids; not facet ids
		void LPUnfolding::find_conflict_loops(const vector<float>& weight, list< list<uint> >& loops)
		{
			int overlap_count = this->m_unfolder->buildFromWeights(weight);
			auto mesh = this->m_unfolder->getModel();
			this->m_unfolder->rebuildModel(); //is this needed?

			vector< vector<uint> > children;
			vector<uint> overlap_faces;
			process_unfolded_tree(children, overlap_faces);

			const auto& collisions = m_unfolder->getOverlppingFacePairs();

			for (uint fid : overlap_faces)
			{
				const auto& cd = collisions[fid];
				for (auto of : cd) //overap face
				{
					if (fid>of) continue; //this will be handled by (of, fid) pair
					//find path connecting fid to of
					list<uint> path;
					find_path(fid, of, children, path);
					if (path.empty() == false)
					{
						//convert the path to edge ids
						list<uint> eids;
						bool r = fids_to_eids(path, eids);
						if(r) loops.push_back(eids);
					}
					else {
						cerr << "! Error: LPUnfolding::find_conflict_loops: Cannot find a path connecting faces: " << fid << " and " << of << endl;
					}
				}//end for of
			}//end for fid
		}


		//find a shortest path connecting f1 and f2
		void LPUnfolding::find_path(uint f1, uint f2, const vector< vector<uint> >& children, list<uint>& path)
		{

			model * m = m_unfolder->getModel();
			//const auto& collisions = m_unfolder->getOverlppingFacePairs();
			std::deque<uint> open;

			//record what has been included in this CC
			//  bool * visited = new bool[m->t_size];
			// memset(visited, false, sizeof(bool) * m->t_size);

			memset(m_visited, 0, sizeof(unsigned char)* m->t_size);

			//add root (f1)
			m_visited[f1] = true;
			open.push_back(f1);
			m->tris[f1].parent_id = f1;

			//loop until all faces in open are handled
			while (open.empty() == false)
			{
				uint fid = open.front();
				open.pop_front();

				if (fid == f2) //path found!
				{
					break; //done
				}

				vector<uint> next = children[fid];
				uint pid = m_unfolder->getParentFaceId(fid);
				if (pid != -1)
					next.push_back(pid);

				for (const auto kid : next)
				{
					if (m_visited[kid]) continue;
					open.push_back(kid);
					m_visited[kid] = true;
					m->tris[kid].parent_id = fid;
				} //end for kid

			} //end while

			int now = f2;
			while (now != f1)
			{
				path.push_front(now);
				now = m->tris[now].parent_id;
			}
			path.push_front(f1);
		}

		/// children points of the unfolded tree
		void LPUnfolding::process_unfolded_tree(vector< vector<uint> >& children, vector<uint>& overlap_faces) const
		{
			//cout << "-----------------" << endl;

			model * m = m_unfolder->getModel();
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

		//convert a path of connected face ids to edge ids
		bool LPUnfolding::fids_to_eids(const list<uint>& fids, list<uint>& eids) const
		{
			if (fids.size() < 2) return false;

			auto ptr = fids.begin();
			auto next = ptr; next++;
			auto mesh = this->m_unfolder->getModel();

			while (next != fids.end())
			{
				//find shared edge between ptr and next
				const auto& tri = mesh->tris[*ptr];
				bool found = false;
				for (short d = 0; d < 3; d++)
				{
					const uint eid = tri.e[d];
					const auto& e = mesh->edges[eid];
					if (e.fid[0] == *next || e.fid[1] == *next)
					{
						found = true;
						eids.push_back(eid);
						break;
					}
				}
				//
				if (found == false)
				{
					cerr << "! Error: LPUnfolding::fids_to_eids: Failed to convert face ids to edge ids. Faces are not connected?" << endl;
					return false;
				}
				//
				ptr = next;
				next++;
			}

			return true;
		}

		uint LPUnfolding::get_neighbor_triangle(model * mesh, uint fid, short side) const
		{
			vector<uint> nei_fids; //neighbor face ids
			triangle& t = mesh->tris[fid];
			edge& e = mesh->edges[t.e[side]];
			//vector<uint> f = e.fid;
			uint nei_fid = (fid == e.fid[0]) ? e.fid[1] : e.fid[0];

			return nei_fid;
		}

	}// end namespace unfolding
}//end namespace masc
