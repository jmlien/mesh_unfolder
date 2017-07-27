#include "polygon/net_analyzer.h"
#include "polygon/chull.h"
#include "polygon/bpc.h"

namespace masc{
	namespace unfolding {
		namespace util {

		void NetAnalyzer::analyze(vector<int>& boundary_vertices, vector<Point2d>& unfolded_positions)
		{

			border_cuts_.clear();//remove existing data

			//convert the boundary to a polygon
			BuildBoundaryPolygon(boundary_vertices, unfolded_positions);

			//build convex hull
			list<polygon::ply_vertex*> hull;
			polygon::hull2d(this->net_.getHead(), this->net_.getHead(), hull);
      for (auto v : hull) v->getExtra().on_hull=true;

			//build bridge
			set<uint> border_edges; //this is a set of edges that are considered hard to fold/glue
      //
			// we find cut edges that do not have its cut halfedges in the same pocket
			//
			for (auto it = hull.begin(); it != hull.end();it++)
			{
				auto nit = it; nit++;
				if (nit == hull.end()) nit = hull.begin();
				if (*it == *nit) continue; //same node?

				polygon::c_BPC bpc;
				if (bpc.build(*it, *nit) == false) continue; //failed for some reason
				auto bridge = bpc.getBridgeEnds();
				//check if the cut edge pairs are in the same pocket
        polygon::ply_vertex* s=bpc.getSource1();
        polygon::ply_vertex* t=bpc.getSource2();

				unordered_map<uint, short> cutedge_counts; //first element is edge id and the second is count
				for (auto i = s; i != t; i = i->getNext())
				{
					cutedge_counts[i->getExtra().mesh_eid]++;
				}//end for i

				//check the count
				for (auto count : cutedge_counts)
				{
					if (count.second == 2) continue; //good edge
					border_edges.insert(count.first);
					//cout << "border cut edge = " << count.first << endl;
				}//end for i

			}//end for it

			//identify zip edges and remove them from "border_edges"
			polygon::ply_vertex * ptr=this->net_.getHead();
			do {
				auto next=ptr->getNext();
				auto pre=ptr->getPre();
				list<uint> matches, vids;
        if(matches.empty()) findMatchingPairs(ptr,next,matches,vids);
        if(matches.empty()) findMatchingPairs(pre,ptr,matches,vids);

        //remove these nice zip vertices/edges from border_edges
				for(uint eid:matches) border_edges.erase(eid);

				//remember the results
				if(matches.empty()==false){
					 this->zip_vids_.push_back(vids);
					 this->zip_eids_.insert(this->zip_eids_.end(),matches.begin(),matches.end());
				 }

				ptr=next;
			} while(ptr!=this->net_.getHead());

			//Expand zip edges
			//TODO: Yun-hyeong will work on this. (remove this line after taks is done)
			// (1) loop around the boundary, for every list of consecutive zip edges, expend the list when possible
			// (2) add the new zip edges to this->zip_eids_
			// (3) repeat (1) and (2) until no new zip edges can be found

			//done, save
			border_cuts_.insert(border_cuts_.end(),border_edges.begin(),border_edges.end());

		}//end NetAnalyzer::analyze

		// Find the boundary for the given model. (copied from SVGWriter.h/cpp)
		void NetAnalyzer::BuildBoundaryPolygon(vector<int>& boundary_vertices, vector<Point2d>& pts)
		{
			auto bv_id = boundary_vertices.begin();

			this->net_.destroy();
			this->net_.beginPoly();
			//cout<<"*bv_id=";
			for (const auto& pt : pts)
			{
				auto v=this->net_.addVertex(pt[0], pt[1]);
				v->getExtra().mesh_vid = *bv_id;
				//cout<<*bv_id<<", ";
				bv_id++;
			}
			this->net_.endPoly();
			//cout<<endl;

			//collect edge ids
			//cout << "eid=";
			auto v = this->net_.getHead();
			do{
				const vertex& meshv = this->model_->vertices[v->getExtra().mesh_vid];
				auto nv = v->getNext();
				uint n_vid = nv->getExtra().mesh_vid;

				for (uint eid : meshv.m_e)
				{
					const edge& e = this->model_->edges[eid];
					if (e.vid[0] == n_vid || e.vid[1] == n_vid)
					{
						v->getExtra().mesh_eid = (e.parent_id == UINT_MAX) ? eid : e.parent_id;
						break;
					}
				}//end for eid

				//cout << v->getExtra().mesh_eid << ", ";
				v = nv;
			}
			while (v != this->net_.getHead());
			//cout << endl;
		}
		// find matching pairs
		void NetAnalyzer::findMatchingPairs
		(polygon::ply_vertex * pre, polygon::ply_vertex * next, list<uint>& matches, list<uint>& vids)
		{
			  auto eid=pre->getExtra().mesh_eid;
				while(eid==next->getExtra().mesh_eid && pre!=next)
				{
					  matches.push_back(eid);

						vids.push_front(pre->getExtra().mesh_vid);
						vids.push_back(next->getExtra().mesh_vid);

						pre=pre->getPre();
						next=next->getNext();
						eid=pre->getExtra().mesh_eid;

						if(eid==matches.front()) break; //visited...
				}//end while

				if(matches.empty()==false) vids.push_back(next->getExtra().mesh_vid);
		}

		}//end namespace util
	}//end namespace unfolding
}//end namespace masc
