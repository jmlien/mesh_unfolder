#pragma once

#include "model.h"
#include "config.h"
#include "polygon/polygon.h"

namespace masc{
	namespace unfolding {
		namespace util {

			class NetAnalyzer
			{
			public:
				NetAnalyzer(const model* model, const Config& config) : model_(model), config_(config), net_(polygon::c_ply::POUT)
				{

				}

				~NetAnalyzer() {

				}

				void analyze(vector<int>& boundary_vertices, vector<Point2d>& unfolded_positions);

        //access Functions

        //get a list of cut edges that are hard to put together
        const list<uint>& getBorderCutEdges() const { return  border_cuts_; }

        //get a list of vertices that can be easily glued (zipped)
        const list< list<uint> >& getZipVertices() const { return zip_vids_; }
        const list<uint>& getZipEdges() const { return zip_eids_; }

        //get net
        const polygon::c_ply & getNetPolygon() const { return net_; }

			protected:

				// builid polygon
				void BuildBoundaryPolygon(vector<int>& boundary_vertices, vector<Point2d>& pts);

        // find matching pairs
        void findMatchingPairs(polygon::ply_vertex * pre, polygon::ply_vertex * next, list<uint>& matches, list<uint>& vids);

			private:

				// Do not own the objective, must be alive during the life-cycle of the SVGWriter object.
				const model* model_;

				// Do not own the objective.
				const Config& config_;

				//a net boundary
				polygon::c_ply net_;

        //identified border cuts
        //these are cuts that from different pockets and cannot be continuesly glued from neighbors
        //these are supposed to be hard to glue vertices/edges
        list<uint> border_cuts_;

        //these are a list of vertices that are easy to glue after cut
        //each list composed of vertices that can be glued in sequence
        //the frist element of a list<uint> is a vertex that has only one cut edge
        //and can be
        list< list<uint> > zip_vids_;
				list<uint> zip_eids_; //this is a list original eids of 3D mesh on zip line
			};

		}//end namespace util
	}//end namespace unfolding
}//end namespace masc
