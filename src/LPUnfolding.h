
/*
* LPUnfolding.h
*
*  Created on: Feb 29, 2015
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

		// Linear Programming Unfolding using branch and abound method
		class LPUnfolding
		{
		public:
		public:

			LPUnfolding(Unfolder* unfolder);
			virtual ~LPUnfolding();


			bool setup(const string& filename);

			//////////////////////////////////////////
			// run n unfoldings and cluster
			//////////////////////////////////////////

			virtual void run();

			//////////////////////////////////////////

			virtual void print(ostream& out) const;


			//////////////////////////////////////////
			// given a weight, unfold the model, find the overlaps. 
			// for each pair of overlapping facets, find a shortest path (i.e loop) connecting them.
			// note that each loop is composed of a list of edge ids; not facet ids
			void find_conflict_loops(const vector<float>& weight, list< list<uint> >& loops);

			// given a weight (1 means cut, 0 means no cut), find the connected components. 
			// note that each loop is composed of a list of edge ids shared between different compnents; 
			// not facet ids
			int find_face_ccs(const vector<float>& weight, list< list<uint> >& ccs);


			//////////////////////////////////////////
			// access
			//////////////////////////////////////////

		protected:

			// will be called once after setup,
			// can be used to override parameters
			void virtual init();

			// setup the entire problem from stream
			virtual bool setup(istream& in);

			// setup the problem from tokens
			virtual bool setup(vector<string>& tokens);

			struct LP_constraints
			{
				LP_constraints(){ type = 1;  upper_bound = lower_bound = 0; }
				vector<uint> eids; //edges involed in this constraint

				int type;  //GLP_FR    free (unbounded) variable, (1) 
			 			   //GLP_LO    lower bound      
						   //GLP_UP    upper bound        
						   //GLP_DB    double bound       
						   //GLP_FX    fixed         

				int upper_bound;
				int lower_bound;
			};


			//sovle the LP problem
			//return true when optimal solution is found
			bool SolveLP(model * m, vector<LP_constraints>& constaints, vector<float>& solution);

			void GenerateConstraints(model * m, vector<LP_constraints>& constaints);
			int AddLazyConstraints(model * m, vector<float>& weight, vector<LP_constraints>& constaints);

			void process_unfolded_tree(vector< vector<uint> >& children, vector<uint>& overlap_faces) const;

			//get the connect component that contains the given facet "start"
			//in weight: 1 mean cut and 0 mean no cut
			void find_face_cc(uint start, const vector<float>& weight, list<uint> & cc);

			//find a path in the unfolding tree
			void find_path(uint f1, uint f2, const vector< vector<uint> >& children, list<uint>& path);

			//convert a path of connected face ids to edge ids
			//return false when failed
			bool fids_to_eids(const list<uint>& fids, list<uint>& eids) const;

			//get the neighbor triangle id...
			//this should really be moved model
			uint get_neighbor_triangle(model * mesh, uint fid, short side) const;

			Unfolder* m_unfolder;
			Unfolder* m_unfolder_from_caller;
			unsigned char* m_visited;
		};
	}
}