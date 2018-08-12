/*
 * UnfoldingEvaluator.cpp
 *
 *  Created on: Feb 29, 2016
 *      Author: zxi
 */

#include "UnfoldingEvaluator.h"

#include <numeric>

#include "libga/Individual.h"
#include "unfolder.h"
#include "OverlappingChecker.h"
#include "util/UnfolderHelper.h"

#include "polygon/polygon.h"
#include "CurveMatching/img2ply.h"
#include "CurveMatching/opencvstd.h"
#include "CurveMatching/curve_db_param.h"
#include "CurveMatching/shadow.h"

namespace masc {
namespace unfolding {

///////////////////////////////////////////////////////////////////////////////
// OverlappingEvaluator
///////////////////////////////////////////////////////////////////////////////
double OverlappingEvaluator::evaluate(masc::ga::Individual* ind) {

  // Build from weights, check global overlaps
  const int global_overlaps = this->m_unfolder->buildFromWeights(ind->getGenome());

  // Check local overlaps
  // local overlaps is not very useful when the net is large.
  const int local_overlaps = this->m_unfolder->checkLocalOverlaps();

  // use default overlapping check
  double fitness = -(global_overlaps + local_overlaps * 100);

  return fitness;
}

// /*
// ///////////////////////////////////////////////////////////////////////////////
// // Learning Evaluator
// ///////////////////////////////////////////////////////////////////////////////
// double LearningEvaluator::evaluate(masc::ga::Individual* ind) {
//
//     // Build from weights, check global overlaps
//     const auto global_overlaps = this->m_unfolder->buildFromWeights(
//                                                                     ind->getGenome());
//
//     // Check local overlaps
//     const auto local_overlaps = this->m_unfolder->checkLocalOverlaps();
//
//     //if (m_params->count("w_overlap")) {
//     //fitness -= (global_overlaps + local_overlaps * 1e6)
//     //* atof(m_params->find("w_overlap")->second.c_str());
//     //}
//
//     // use default overlapping check
//     double fitness = -(global_overlaps + local_overlaps * 100);
//
//     const auto& config = this->m_unfolder->getConfig();
//
//     this->m_unfolder->rebuildModel();
//
//     if (this->m_unfolder->isFlattened()) {
//
//         if (m_params->count("w_bias"))
//         fitness = 100;//atof(m_params->find("w_bias")->second.c_str());
//         else
//         fitness = 100;
//         int w_count = 0;
//
//         util::UnfolderHelper unfolder_helper(m_unfolder);
//
//         int num_leaf_nodes=0;
//         int num_1_degree_nodes=0;
//         int num_2_degree_nodes=0;
//         int num_3_degree_nodes=0;
//         int total_degrees=0;
//
//         for (int i = 0; i < unfolder_helper.m_num_children.size(); i++) {
//
//             switch(unfolder_helper.m_num_children[i]) {
//                 case 0: num_leaf_nodes++; break;
//                 case 1: num_1_degree_nodes++; break;
//                 case 2: num_2_degree_nodes++; break;
//                 case 3: num_3_degree_nodes++; break;
//             }
//             total_degrees += unfolder_helper.m_num_children[i];
//
//         }
//
//
//         if (m_params->count("w_cutlength")) {
//             //cout << "cutlength = " << unfolder_helper.getTotalCutLengthNormalized() << endl;
//             fitness += unfolder_helper.getTotalCutLengthNormalized()
//             * atof(m_params->find("w_cutlength")->second.c_str());
//             w_count++;
//         }
//
//         if (m_params->count("w_hullarea")) {
//             //cout << "huallarea = " << this->m_unfolder->getHullArea() << endl;
//
//             fitness += this->m_unfolder->getModel()->surface_area
//             / this->m_unfolder->getHullArea()
//             * atof(m_params->find("w_hullarea")->second.c_str());
//             w_count++;
//         }
//
//
//         if (m_params->count("w_num_leaf_nodes")) {
//             //cout << "num_leaf_nodes = " <<  (double) num_leaf_nodes / (double)total_degrees << endl;
//             fitness += (double) num_leaf_nodes / (double)total_degrees
//             * atof(m_params->find("w_num_leaf_nodes")->second.c_str());
//             w_count++;
//         }
//
//         if (m_params->count("w_num_1_degree_nodes")) {
//             //cout << "num_1_degree_nodes = " << (double) num_1_degree_nodes / (double)total_degrees << endl;
//             fitness += (double) num_1_degree_nodes / (double)total_degrees
//             * atof(m_params->find("w_num_1_degree_nodes")->second.c_str());
//             w_count++;
//         }
//
//         if (m_params->count("w_num_2_degree_nodes")) {
//             //cout << "num_2_degree_nodes = " << (double) num_2_degree_nodes / (double)total_degrees << endl;
//             fitness += (double) num_2_degree_nodes / (double)total_degrees
//             * atof(m_params->find("w_num_2_degree_nodes")->second.c_str());
//             w_count++;
//         }
//
//         if (m_params->count("w_border_cuts")) {
//             //cout << " border_cuts " << (double)unfolder_helper.computeBorderCutsLength() << endl;
//             fitness += (double)unfolder_helper.computeBorderCutsLength()
//             / (double) unfolder_helper.getTotalCutLength()
//             * atof(m_params->find("w_border_cuts")->second.c_str());
//             w_count++;
//         }
//         //cout << "Learning " << w_count << " parameters ..." << endl;
//         //cout << "Learning fittness = " << fitness << endl;
//         //cout << "n1d = " << num_leaf_nodes << endl;
//         //cout << "n2d = " << num_1_degree_nodes << endl;
//
//     }
//
//     return fitness;
// }*/


double LearningEvaluator::evaluate(masc::ga::Individual* ind) {

    // Build from weights, check global overlaps
    const auto global_overlaps = this->m_unfolder->buildFromWeights(
                                                                        ind->getGenome());
    // Check local overlaps
    const auto local_overlaps = this->m_unfolder->checkLocalOverlaps();

    //if (m_params->count("w_overlap")) {
    //fitness -= (global_overlaps + local_overlaps * 1e6)
    //* atof(m_params->find("w_overlap")->second.c_str());
    //}

    // use default overlapping check
    double fitness = -(global_overlaps + local_overlaps * 100);
    double diff = 0.0;

    const auto& config = this->m_unfolder->getConfig();

    this->m_unfolder->rebuildModel();

    if (this->m_unfolder->isFlattened()) {

        if (m_params->count("w_bias"))
            fitness = atof(m_params->find("w_bias")->second.c_str());
        else
            fitness = 0.0;
        int w_count = 0;

        double SCALE = 1.0;

        util::UnfolderHelper unfolder_helper(m_unfolder);

        int num_leaf_nodes=0;
        int num_1_degree_nodes=0;
        int num_2_degree_nodes=0;
        int num_3_degree_nodes=0;
        int total_degrees=0;

        for (int i = 0; i < unfolder_helper.m_num_children.size(); i++) {

            switch(unfolder_helper.m_num_children[i]) {
            case 0: num_leaf_nodes++; break;
            case 1: num_1_degree_nodes++; break;
            case 2: num_2_degree_nodes++; break;
            case 3: num_3_degree_nodes++; break;
            }
            total_degrees += unfolder_helper.m_num_children[i];

        }

        if (m_params->count("w_cutlength")) {
            //cout << "cutlength = " << unfolder_helper.getTotalCutLengthNormalized() << endl;
            diff = unfolder_helper.getTotalCutLengthNormalized()
                - atof(m_params->find("w_cutlength")->second.c_str());
            fitness = fitness-diff*diff*SCALE+2.0*SCALE;
            w_count++;
        }

        if (m_params->count("w_hullarea")) {
            //cout << "huallarea = " << this->m_unfolder->getModel()->surface_area / this->m_unfolder->getHullArea() << endl;

            diff = this->m_unfolder->getModel()->surface_area
                / this->m_unfolder->getHullArea()
                - atof(m_params->find("w_hullarea")->second.c_str());
            fitness = fitness-diff*diff*SCALE+2.0*SCALE;
            w_count++;
        }

        if (m_params->count("w_num_leaf_nodes")) {
            //cout << "num_leaf_nodes = " <<  (double) num_leaf_nodes / (double)total_degrees << endl;
            diff = (double) num_leaf_nodes / (double)total_degrees
                - atof(m_params->find("w_num_leaf_nodes")->second.c_str());
            fitness = fitness-diff*diff*SCALE+2.0*SCALE;
            w_count++;
        }

        if (m_params->count("w_num_1_degree_nodes")) {
            //cout << "num_1_degree_nodes = " << (double) num_1_degree_nodes / (double)total_degrees << endl;
            diff = (double) num_1_degree_nodes / (double)total_degrees
                - atof(m_params->find("w_num_1_degree_nodes")->second.c_str());
            fitness = fitness-diff*diff*SCALE+2.0*SCALE;
            w_count++;
        }

        if (m_params->count("w_num_2_degree_nodes")) {
            //cout << "num_2_degree_nodes = " << (double) num_2_degree_nodes / (double)total_degrees << endl;
            diff = (double) num_2_degree_nodes / (double)total_degrees
                - atof(m_params->find("w_num_2_degree_nodes")->second.c_str());
            fitness = fitness-diff*diff*SCALE+2.0*SCALE;
            w_count++;
        }

        if (m_params->count("w_border_cuts")) {
            //cout << " border_cuts " << (double)unfolder_helper.computeBorderCutsLength() / (double) unfolder_helper.getTotalCutLength() << endl;
            diff = (double)unfolder_helper.computeBorderCutsLength()
                / (double) unfolder_helper.getTotalCutLength()
                - atof(m_params->find("w_border_cuts")->second.c_str());
            fitness = fitness-diff*diff*SCALE+2.0*SCALE;
            w_count++;
        }
        //cout << "Learning " << w_count << " parameters ..." << endl;
        //cout << "Learning fittness = " << fitness << endl;
        //cout << "n1d = " << num_leaf_nodes << endl;
        //cout << "n2d = " << num_1_degree_nodes << endl;

    }
    return fitness;
}

///////////////////////////////////////////////////////////////////////////////
// AreaEvaluator
///////////////////////////////////////////////////////////////////////////////
AreaEvaluator::AreaEvaluator(Unfolder* unfolder) :
    UnfoldingEvaluator(unfolder) {
  m_checker.reset(new PixelChecker(unfolder));
  m_best_ratio = 1e3;
}

double AreaEvaluator::evaluate(masc::ga::Individual* ind) {
  // build the model, but do not check overlaps
  this->m_unfolder->buildFromWeights(ind->getGenome(), false);

  // check area ratio
  const double area_ratio = this->m_checker->checkOverlapping(
      m_unfolder->getUnfolded(), m_unfolder->getConfig());

  // maximum overlaps
  int overlaps = (m_unfolder->getModel()->t_size)
      * (m_unfolder->getModel()->t_size);

  // only check when ratio is low
  if (area_ratio < m_best_ratio * 1.01)
    overlaps = m_unfolder->checkOverlaps();

  if (area_ratio < m_best_ratio)
    m_best_ratio = area_ratio;

  return -overlaps;
}


///////////////////////////////////////////////////////////////////////////////
// Net Evaluators
///////////////////////////////////////////////////////////////////////////////

double CutLengthEvaluator::evaluate(Unfolder* unfolder) {
  return unfolder->getTotalCutLength();
}

double HullAreaEvaluator::evaluate(Unfolder* unfolder) {
  return 1.0/unfolder->getHullArea();
}

// estimate max side length of bounding box
double BoxSideEvaluator::evaluate(Unfolder* unfolder) {

  float min_x=FLT_MAX, max_x=-FLT_MAX, min_y=FLT_MAX, max_y=-FLT_MAX;

  const MESH& unfolded_mesh=unfolder->getUnfolded();
  for (const auto& f : unfolded_mesh) {
    for (int i = 0; i < 3; ++i)
    {
      const Vector3d& pos=f[i].second;
      if(pos[0]<min_x) min_x=pos[0];
      if(pos[0]>max_x) max_x=pos[0];
      if(pos[2]<min_y) min_y=pos[2];
      if(pos[2]>max_y) max_y=pos[2];
    }
  }

  return unfolder->getTotalCutLength()/max((max_x-min_x), (max_y-min_y));
}

//-----------------------------------------------------------------------------
//
// parameters for PolygonFitEvaluator
//
//-----------------------------------------------------------------------------

//a shape is said to be inside the target is the boundary and area differences are smaller than the values below
const int   matching_tolerable_max_diff = 3;   //max boundary difference
const float matching_tolerable_sum_diff = 200; //max area difference

//smallest boundary that will be extracted from the image
const int curveDB_resample_size=100;
const int spot_target_smallest_curve_size = 70; // 5, 10, 25
const int spot_target_longest_curve_size = 99;
const int spot_target_offset_size = 2;

typedef CurveMatcher MATCHER;
typedef cv::Point CVPT;
typedef CURVE_DB_PARAM<MATCHER> MY_CURVE_DB_PARAM;

PolygonFitEvaluator::PolygonFitEvaluator(const std::string& stencil_filename)
:NetEvaluator(), m_poly_ptr(NULL), m_best_net_ptr(NULL), m_min_error(FLT_MAX)
{
  auto * target = new CSShape<MATCHER, CVPT>();
  assert(target);

  if(stencil_filename.empty())
  {
    cerr<<"! Error: No stencil file is given for PolygonFitEvaluator"<<endl;
    exit(1);
  }

  cv::Mat img = cv::imread(stencil_filename, CV_LOAD_IMAGE_GRAYSCALE);

  //cout<<"image size="<<img.size[0]<<" "<<img.size[1]<<endl;

  // flip the image
  //cv::flip(img, img, 0);

  //cv::imshow("img",img);
  //waitKey(0);
  //
  c_polygon polygon;
  bool success = img2ply(img, polygon);
  if(!success)
  {
    cerr<<"! Error: Failed to create polygon from stencil file "<<stencil_filename<<" in PolygonFitEvaluator"<<endl;
    exit(1);
  }

  //cout<<polygon<<endl;
  target->contours.push_back(polygon);

  //from polygon create matching curves
  MY_CURVE_DB_PARAM::CurveResampleSize = curveDB_resample_size;
  MY_CURVE_DB_PARAM::SmallestCurveSegmentSize = spot_target_smallest_curve_size;
  MY_CURVE_DB_PARAM::LongestCurveSegmentSize = spot_target_longest_curve_size;
  MY_CURVE_DB_PARAM::OffsetStepSize = spot_target_offset_size;

  target->buildCurveSegmentDB();
  //target->image = img;//cv::Mat(img.size[0], img.size[1], CV_8UC1, img);
  //m_target.build_domain(shadow_shader, depthVP, buffers);
  //target->builDistField();
  m_poly_ptr=target;
}


PolygonFitEvaluator::~PolygonFitEvaluator()
{
  auto * target = (CSShape<MATCHER, CVPT>*)m_poly_ptr;
  auto * source = (CSShape<MATCHER, CVPT>*)m_best_net_ptr;

  if(target!=NULL && source!=NULL)
  {
    auto& target_db = target->contourDBs[0]; //only consider the first contour
    auto& source_db = source->contourDBs[0]; //only consider the first contour

#if 1
      MATCHER matcher;
      MY_CURVE_DB_PARAM::initializeMatcher(matcher);
      matcher.setTarget(target_db.cv_contour);
      matcher.setSource(source_db.cv_contour);

      //save the best match to file
      //cout<<"min_error="<<min_error<<endl;
      //matcher.visualizeMatching(best_source_curve_segment.x,best_source_curve_segment.y, best_target_curve_segment.x, best_target_curve_segment.y);
      stringstream ss;
      ss<<"polygonfitevaluator_best_net_gen_err_"<<m_min_error<<".jpg";
      matcher.renderMatching(ss.str(), best_src_off,best_src_len, best_target_off, best_src_len);
      cout<<"- PolygonFitEvaluator:: Saved best matching to "<<ss.str()<<endl;
#endif
  }

  delete target;
}

//convert a flattened model into a CSShape
inline void model2poly(const model* m, CSShape<MATCHER, CVPT>& net, const Config& config)
{
  SVGWriter writer(m, config);
  vector<int> boundary;
  writer.FindBoundaryPolygon(&boundary);

  //create netshadow
  c_polygon polygon;
  c_ply ply(c_ply::POUT);
  ply.beginPoly();
  for(auto& id : boundary)
  {
    Vector3d pos=writer.GetSVGCoord(id);
    ply.addVertex(pos[0],pos[2]);
  }
  ply.endPoly();
  polygon.push_back(ply);
  net.contours.push_back(polygon);
  net.buildCurveSegmentDB();
}

double PolygonFitEvaluator::evaluate(Unfolder* unfolder)
{
  //cout<<"is flattended="<<unfolder->isFlattened()<<endl;

  CSShape<MATCHER, CVPT> net;
  auto * target = (CSShape<MATCHER, CVPT>*)m_poly_ptr;
  if(target==NULL)
  {
      cerr<<"! Error: Cannot convert from m_poly_ptr to CSShape"<<endl;
      exit(1);
  }

  //auto start_time = clock();

  unfolder->rebuildModel();
  model2poly(unfolder->getNet(), net, unfolder->getConfig());

  //cout<<"model2poly takes "<<(clock()-start_time)*1.0/CLOCKS_PER_SEC<<" secs"<<endl;

  //compute transforms and check if we can put the net inside...
  auto& target_db = target->contourDBs[0]; //only consider the first contour
  auto& source_db = net.contourDBs[0]; //only consider the first contour
  int cs_size = target_db.curve_segments.size();

  //initialize the matcher
  MATCHER matcher;
  MY_CURVE_DB_PARAM::initializeMatcher(matcher);
  matcher.setTarget(target_db.cv_contour);

  //initalize the matches
  //update the best match so far
  auto best_target_curve_segment = target_db.curve_segments[0];
  auto best_source_curve_segment = source_db.curve_segments[0];
  double min_error=FLT_MAX;

  //cout<<"target cs_size="<<cs_size<<" source cs_size="<<source_db.curve_segments.size()<<endl;

  //start_time = clock();

  for (int j = 0; j < cs_size; j++)
  {

    //cout << "\t- Working on curve segetment=" << j << "/" << cs_size << " updated!!" << endl;
    auto& target_curve_segment = target_db.curve_segments[j];
    vector<double>& target_curve_curvatures = target_db.curvatures[j];

    matcher.setSource(source_db.cv_contour);

    cv::DMatch match;

    //match the shadow to the segment of the target...
    //the return is the best match
    //matcher.CompareCurvesUsingSignatureDB(source_db.curve_segments, target_curve_segment, source_db.curvatures, target_curve_curvatures, match);
    matcher.CompareCurvesUsingSignatureDB_curvature_only(source_db.curve_segments, target_curve_segment, source_db.curvatures, target_curve_curvatures, match);

    if (match.queryIdx < 0) continue; //there is not match to the curve_segment
    //

    //make sure that the entire contour can fit in....
		// int shadow_len = source_db.curve_segments[match.queryIdx].x;
		// int shadow_off = source_db.curve_segments[match.queryIdx].y;
		// int target_len = target_curve_segment.x;
		// int target_off = target_curve_segment.y;

    //double rmse = matcher.ComputeWholeRMSE(shadow_len, shadow_off, target_len, target_off);
    //if (rmse >= min_error) continue; //the error is bigger
    if(match.distance>=min_error) continue; //the error is bigger

    //cv::Mat transform;
    //matcher.computeTransform(shadow_len, shadow_off, target_len, target_off, transform);
    //bool isinside = target->inside(source_db.cv_contour, transform, matching_tolerable_max_diff, matching_tolerable_sum_diff);

    //if (isinside == false) continue; //does not fit to the target...

    best_target_curve_segment = target_db.curve_segments[j];
    best_source_curve_segment = source_db.curve_segments[match.queryIdx];
		//min_error = rmse;
    min_error = match.distance;
  }

   //cout<<"find match takes "<<(clock()-start_time)*1.0/CLOCKS_PER_SEC<<" secs"<<endl;
   if(min_error<m_min_error)
   {
      m_min_error=min_error;
      best_src_off=best_source_curve_segment.x;
      best_src_len=best_source_curve_segment.y;
      best_target_off=best_target_curve_segment.x;
      best_target_len=best_target_curve_segment.y;
      if(m_best_net_ptr==NULL)  m_best_net_ptr= new CSShape<MATCHER, CVPT>();
      assert(m_best_net_ptr);
      //copy net to m_best_net_ptr
      *((CSShape<MATCHER, CVPT>*)m_best_net_ptr)=net;
   }

#if 0
  //save the best match to file
  //cout<<"min_error="<<min_error<<endl;
  //matcher.visualizeMatching(best_source_curve_segment.x,best_source_curve_segment.y, best_target_curve_segment.x, best_target_curve_segment.y);
  static int counter=0;
  stringstream ss;
  ss<<"gen_err_"<<min_error<<"_"<<counter++<<".jpg";
  matcher.renderMatching(ss.str(), best_source_curve_segment.x,best_source_curve_segment.y, best_target_curve_segment.x, best_target_curve_segment.y);
  //cout<<"\t- Saved best matching to "<<ss.str()<<endl;
#endif

  return 1.0/min_error;
}


} /* namespace unfoldings */
} /* namespace masc */
