/*
 * config.h
 *
 *  Created on: Nov 21, 2014
 *      Author: zxi
 */

#ifndef CONFIG_H_
#define CONFIG_H_

#include <ctime>
#include <climits>

#include <string>
using namespace std;

#include "mathtool/Vector.h"
using namespace mathtool;

enum class CutHeuristic {
  FLAT_TREE,          // f
  UNFLAT_TREE,        // uf
  MINIMUM_PERIMETER,  // p
  MAXIMUM_PERIMETER,  // mp
  STEEPEST_EDGE,      // s
  UNSTEEPEST_EDGE,    // us
  RANDOM,             // r
  RANDOM_EDGE,        // re
  BRUTE_FORCE,        // b
  GA,                 // ga
  CLUSTERING,         // c
  LP,                 // lp
  BLOOMING,           // bloom
  LASER,              // laser
  EASILY_FOLDABLE,    // easy
};

enum class Objective {
  CUT_LENGTH,         // total cut length
  HULL_AREA,          // hull area
  BOX_SIDE,           // max side length of bounding box
  POLYGON_FIT,        // polygon fit
};

//net stitching and cutting
enum class NetSurgery
{
  NO_SURGERY,             //disalbe net surgery
  SET_COVER_SURGERY,      //split overlapping part of the nete into disjoint parts
  TOPOLOGICAL_SURGERY,    //extension of set cover that merges split nets using GA
  SUBDIVID_SURGERY,       //subdivid a triangle and rearrange sub-triangles to avoid overlap
  BOXING_SURGERY,         //split a net into one or multiple boxes of given size
  CAGING_SURGERY          //ensure the net is contained in a polygonal container
};

struct Config {
public:
  Config()
  {
    // flat tree is good for non-convex shapes
    this->heuristic = CutHeuristic::FLAT_TREE;

    this->max_retries = 100;
    this->early_stop = false;

    this->random_baseface = false;
    this->less_cuts = false;
    this->seed = time(nullptr);

    this->fold_to_compact_state = false;
    this->find_compact_folding = false;

    // for ga
    this->ga_config_filename = "unfolding.ga";
    this->ga_max_gen = INT_MAX;
    this->ga_local_adj = 0;

  	//for cluster
    this->cluster_config_filename = "unfolding.cluster";
    this->lp_config_filename = "unfolding.lp";

  	//for blooming
  	this->blooming_strategy = "flower";
  	this->blooming_range = 1e-5; //small value to avoid 90 degree walls disconnecting the faces
  	this->blooming_dir = Vector3d(0, 0, 0);

	   //
    this->optimize_unfolding = false;

    this->k = -1;
    this->run = -1;
    this->scalar_score = false;

    this->quite = false;

    this->shrink = true;
    this->find_boundary = true;
    this->shrink_factor = 0.999;

    this->disable_gui = false;

    this->scale = 1.0;

    this->svg_dump_labels = true;

    this->svg_edge_hints = true;

    this->svg_add_tabs = false;

    this->svg_valley_chamfer = false;

    this->svg_valley_chamfer_width_ratio = 1.0;

    this->svg_boxing_width=13;
    this->svg_boxing_height=13;

    this->no_tick = false;

    this->use_rapid = false;

    this->find_best_base_face = true;

    this->unfolding_motion = Linear_Unfolding;

    this->use_user_vector = false;

    this->record_overlap = true;

    this->dump_wrl = false;

    this->label_font_scale = 1.0;

    this->max_iterations = INT_MAX;

    this->pixel_checker = false;

    this->weighted_dist = false;

    this->force_mode = false;

    this->binary_format = false;

    this->score_only = false;

    this->no_dump = false;
    this->no_rebuild = false;

    this->dump_models_to_wrl = false;

    this->output_vertex_score = false;

    this->repair_mode = false;

    this->extra_cuts_x = 0.7;
    this->extra_cuts_y = 0.1;

    this->test_mode = false;
    this->training_mode = false;

    this->grid_height = 4;
    this->grid_width = 4;
    this->piles = 4;
    this->output_file = "";
    this->params_file_path = "";

    this->opt_obj = Objective::CUT_LENGTH;
    this->surgery_method = NetSurgery::NO_SURGERY;
  }

  CutHeuristic heuristic;
  string filename;

  //////////////////////////////////////////////
  //  for GA
  //////////////////////////////////////////////
  string ga_config_filename;
  // maximum generation to override
  int ga_max_gen;
  // maximum overlap face to start using local adjustment.
  int ga_local_adj;

  //////////////////////////////////////////////
  //  for LP
  //////////////////////////////////////////////
  string lp_config_filename;

  ///////////////////////////////////////////////
  // for compactness
  //////////////////////////////////////////////
  bool fold_to_compact_state;
  bool find_compact_folding;
  string compact_config_filename;

  //////////////////////////////////////////////
  //  for clustering
  //////////////////////////////////////////////
  string cluster_config_filename;
  // number of components (k<=0, use value from config)
  int k;
  // number of runs (not used yet)
  int run;
  // maximum iterations
  int max_iterations;
  // use external labels (exp: first level only)
  string label_filename;
  // weighted distance between faces
  bool weighted_dist;
  // export the weight in binary format
  bool binary_format;
  // only do the scoring
  bool score_only;
  // use a scalar as the score for each face
  bool scalar_score;

  string score_filename;
  bool output_vertex_score;
  string vertex_score_filename;

  // in repair mode, load the score and labels, and repair the clusters
  bool repair_mode;

  //////////////////////////////////////////////
  //  for blooming
  //////////////////////////////////////////////
  string blooming_strategy;
  float blooming_range;
  Vector3d blooming_dir;

  //////////////////////////////////////////////
  //  for shadow container
  //////////////////////////////////////////////
  string stencil_filename; //stencil image, a polygon will be extracted from this image

  //////////////////////////////////////////////
  string weights_filename;
  string cuts_filename;

  //////////////////////////////////////////////
  // reduce the cut length
  bool optimize_unfolding;

  // max times to try for a heuristic method with randomness
  int max_retries;

  // stop when found a net
  bool early_stop;

  // whether to use random base face
  int baseface = -1;
  bool random_baseface;

  // less cut on tessellated face
  bool less_cuts;
  bool quite;
  bool shrink;
  bool find_boundary;
  double shrink_factor;

  // scale of the net
  double scale;

  // disable gui, find and dump net
  bool disable_gui;

  // seed to use
  int seed;

  // dump labels in the svg file
  bool svg_dump_labels;

  // dump edge hints to make foldnig easier
  // only dumpped into cut svg file
  bool svg_edge_hints;

  // whether to add tabs or not
  bool svg_add_tabs;

  // whether to add chamfer or not
  // if so, the width of the chamger is related to folding angle
  bool svg_valley_chamfer;
  float svg_valley_chamfer_width_ratio;

  //desired svg doc width and height
  float svg_boxing_width;
  float svg_boxing_height;

  // label font scale
  double label_font_scale;

  // do not add tick in the dumpped filename
  bool no_tick;

  // whether to use rapid for overlapping/collision detection
  bool use_rapid;

  // whether to find the best base face
  bool find_best_base_face;

  enum Unfolding_Motion {Linear_Unfolding, Ordered_Unfolding, Laser_Unfolding, Uniform_Unfolding};

  // whether the unfolding is ordered
  Unfolding_Motion unfolding_motion;

  //
  bool use_user_vector; //if this is true, user vector will be used. otherwise random vector will be used. default: false
  Vector3d user_vector; //user specified vector

  //
  bool record_overlap; /// set true to turn on recording of overlapping pairs.

  // dump
  bool no_rebuild; //
  bool no_dump; // do not dump anything...
  bool dump_wrl;
  bool dump_models_to_wrl;

  // extra cuts
  double extra_cuts_x;
  double extra_cuts_y;

  // using pixel checker to estimate overlapping
  bool pixel_checker;

  // ignore all cache files
  bool force_mode;

  // path for texture
  string texture_path;

  // for optimization
  Objective opt_obj;

  // for net surgery
  NetSurgery surgery_method;

  // test only
  bool test_mode;

  // training mode, no dump svgs only
  bool training_mode;

  // for hamiltonian paths in grid patterns
  int grid_width; // -gw
  int grid_height; // -gh
  int piles; // -piles
  int rows { -1 };
  int cols { -1 };

  string output_file;

  string params_file_path;

};

#endif /* CONFIG_H_ */
