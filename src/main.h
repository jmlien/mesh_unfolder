//------------------------------------------------------------------------------
//  Copyright 2007-2008 by Jyh-Ming Lien and George Mason University
//  See the file "LICENSE" for more information
//------------------------------------------------------------------------------

#ifndef _BF_MINKOWSKI_SUM_H_
#define _BF_MINKOWSKI_SUM_H_

//
//
//
// This is the main header file is shared by m+3d.cpp (single thread) and pm+3d.cpp
// (multi-thread). Functions defined in this file are mainly for text parsing and
// OpenGL rendering
//
//
//
#include <cfloat>

#include <list>
#include <string>
using namespace std;

#include "objReader.h"
#include "model.h"
#include "unfolder.h"
#include "config.h"

#include "ClusterUnfolding.h"
#include "LPUnfolding.h"
#include "UnfoldingProblem.h"
#include "CompactProblem.h"
using namespace masc::unfolding;

#include "util/DataHelper.h"

//-----------------------------------------------------------------------------
// INPUTS
vector<string> filenames;
Config config;

//-----------------------------------------------------------------------------
// Intermediate data
vector<Unfolder*> unfolders;

double R = 0;       //radius
Point3d COM;     //center of mass

//-----------------------------------------------------------------------------
//for random rotation
double rot_theta;
Vector3d rot_vec;
Quaternion current_rot;

//-----------------------------------------------------------------------------
//read M+ from file
bool readfromfile();
void computeCOM_R();

//-----------------------------------------------------------------------------
//
//
//
//  Open GL stuff below
//
//
//-----------------------------------------------------------------------------

#include "draw.h"

//-----------------------------------------------------------------------------
bool parseArg(int argc, char ** argv) {
  if (argc == 1)
    return false;

  for (int i = 1; i < argc; i++) {
    auto const arg = string(argv[i]);
    if (arg == "-r") {
      config.max_retries = std::stoi(argv[++i]);
    } else if (arg == "-test") {
      config.test_mode = true;
    } else if (arg == "-train") {
        config.training_mode = true;
    } else if (arg == "-h") {

      auto const heuristic = string(argv[++i]);
      if (heuristic[0] == 's')
        config.heuristic = CutHeuristic::STEEPEST_EDGE;
      else if (heuristic == "us")
        config.heuristic = CutHeuristic::UNSTEEPEST_EDGE;
      else if (heuristic[0] == 'f')
        config.heuristic = CutHeuristic::FLAT_TREE;
      else if (heuristic == "uf")
        config.heuristic = CutHeuristic::UNFLAT_TREE;
      else if (heuristic[0] == 'p')
        config.heuristic = CutHeuristic::MINIMUM_PERIMETER;
      else if (heuristic == "mp")
        config.heuristic = CutHeuristic::MAXIMUM_PERIMETER;
      else if (heuristic == "r")
        config.heuristic = CutHeuristic::RANDOM;
      else if (heuristic == "re")
        config.heuristic = CutHeuristic::RANDOM_EDGE;
      else if (heuristic[0] == 'b')
        config.heuristic = CutHeuristic::BRUTE_FORCE;
      else if (heuristic == "ga")
        config.heuristic = CutHeuristic::GA;
      else if (heuristic == "c")
        config.heuristic = CutHeuristic::CLUSTERING;
      else if (heuristic == "lp")
        config.heuristic = CutHeuristic::LP;
      else {
        cerr << "!Error! Unknown heuristic type = " << heuristic << endl;
        return false;
      }
    } else if (arg == "-objective") {
      const auto obj = string(argv[++i]);
      if (obj == "cut_length") {
        config.opt_obj = Objective::CUT_LENGTH;
      } else if (obj == "hull_area") {
        config.opt_obj = Objective::HULL_AREA;
      } else if (obj == "polygon") {
        config.opt_obj = Objective::POLYGON_FIT;
      } else {
        cerr << "!Error! Unknown objective type = " << obj << endl;
        return false;
      }
    } else if (arg == "-g") {
      config.disable_gui = true;
    } else if (arg == "-rb") {
      config.random_baseface = true;
    } else if (arg == "-bf" || arg == "--base-face") {
      config.baseface = std::stoi(argv[++i]);
    } else if (arg == "-mg" || arg == "--max-gen") {
      config.ga_max_gen = std::stoi(argv[++i]);
    } else if (arg == "-la" || arg == "--local-adjust") {
      config.ga_local_adj = std::stoi(argv[++i]);
    } else if (arg == "-lc") {
      config.less_cuts = true;
    } else if (arg == "-s") {
      config.seed = std::stoi(argv[++i]);
    } else if (arg == "-q") {
      config.quite = true;
    } else if (arg == "-ns") {
      config.shrink = false;
    } else if (arg == "-nfb") {
      config.find_boundary = false;
    } else if (arg == "-ga") {
      config.ga_config_filename = argv[++i];
    } else if (arg == "-lp") {
      config.lp_config_filename = argv[++i];
    } else if (arg == "-p") {
      config.find_compact_folding = true;
    } else if (arg == "--compact") {
      config.fold_to_compact_state = true;
      config.compact_config_filename = argv[++i];
    } else if (arg == "-c") {
      config.cluster_config_filename = argv[++i];
    } else if (arg == "-optimize") {
      config.optimize_unfolding = true;
    } else if (arg == "-k") {
      config.k = std::stoi(argv[++i]);
    } else if (arg == "-i") {
      config.max_iterations = std::stoi(argv[++i]);
    } else if (arg == "-run") {
      config.run = std::stoi(argv[++i]);
    } else if (arg == "-weights") {
      config.weights_filename = argv[++i];
    } else if (arg == "-cuts") {
      config.cuts_filename = argv[++i];
    } else if (arg == "-sf") {
      config.shrink_factor = stof(argv[++i]);
    } else if (arg == "-scale") {
      config.scale = stof(argv[++i]);
    } else if (arg == "-no_dump") {
      config.no_dump = true;
    } else if (arg == "--no-rebuild") {
      config.no_rebuild = true;
    } else if (arg == "-nl") {
      config.svg_dump_labels = false;
    } else if (arg == "-no-hints") {
      config.svg_edge_hints = false;
    } else if (arg == "-lfs") {
      config.label_font_scale = stof(argv[++i]);
    } else if (arg == "-nt") {
      config.no_tick = true;
    } else if (arg == "-rapid") {
      config.use_rapid = true;
    } else if (arg == "-tab") {
      config.add_tabs = true;
    } else if (arg == "-nbb") {
      config.find_best_base_face = false;
    } else if (arg == "-ordered") {
      config.ordered_unfolding = true;
    } else if (arg == "-pc") {
      config.pixel_checker = true;
    } else if (arg == "-wrl") {
      config.dump_models_to_wrl = true;
    } else if (arg == "-w") {
      config.weighted_dist = true;
    } else if (arg == "-l") {
      config.label_filename = argv[++i];
    } else if (arg == "-f") {
      config.force_mode = true;
    } else if (arg == "-b") {
      config.binary_format = true;
    } else if (arg == "-so") {
      config.score_only = true;
      cout << " - scale_only = true" << endl;
    } else if (arg == "-ss") {
      config.scalar_score = true;
      cout << " - scale_score = true" << endl;
    } else if (arg == "--score") {
      config.scalar_score = true;
      config.score_filename = argv[++i];
    } else if (arg == "--repair") {
      config.repair_mode = true;
    } else if (arg == "-ecx") {
      config.extra_cuts_x = std::stof(argv[++i]);
    } else if (arg == "-ecy") {
      config.extra_cuts_y = std::stof(argv[++i]);
    } else if (arg == "-texture") {
      config.texture_path = argv[++i];
    } else if (arg == "-early_stop") {
      config.early_stop = true;
    } else if (arg == "-ovs") {
      if (config.score_filename.empty()) {
        cout << "score filename can not be empty!" << endl;
      } else {
        config.output_vertex_score = true;
        config.vertex_score_filename = config.score_filename + ".v";
      }
    } else if (arg == "-gw") {
      config.grid_width = std::stoi(argv[++i]);
    } else if (arg == "-gh") {
      config.grid_height = std::stoi(argv[++i]);
    } else if (arg == "-piles") {
      config.piles = std::stoi(argv[++i]);
    } else if (arg == "-paramsfile") {
      config.params_file_path = argv[++i];
    } else if (arg == "-outputfile") {
      config.output_file = argv[++i];
    } else if (arg == "-polygon") {
      config.stencil_filename = argv[++i];
    } else if (arg[0] == '-') {
      cerr << "!Error! Unknown arg = " << arg << endl;
      return false;
    } else {
      filenames.push_back(string(argv[i]));
    }
  }

  return true;
}

void printUsage(char * name) {

  Config default_config;

  //int offset=20;
  cerr << "Usage: " << name << " [options] *.obj/*.off \n\n";

  cerr << "Heuristic Methods\n";
  cerr << "  -h heuristic | use heuristic method\n";
  cerr << "      s        | STEEPEST_EDGE\n";
  cerr << "      f        | FLAT_TREE (default)\n";
  cerr << "     uf        | UNFLAT_TREE\n";
  cerr << "      p        | MINIMUM_PERIMETER\n";
  cerr << "     mp        | MAXIMUM_PERIMETER\n";
  cerr << "      r        | RANDOM\n";
  cerr << "     ga        | Genetic Algorithm\n";
  cerr << "      c        | Clustering based method\n";
  cerr << "     lp        | Linear Programming branch and bound based method\n";
  cerr << "\n";

  cerr << "Clustering\n";
  cerr << "  -c filename  | specify the clustering config file. defualt = '"
      << default_config.cluster_config_filename << "'\n";
  cerr << "  -k           | override the number of components.\n";
  cerr << "  -i           | maximum iterations.\n";
  cerr << "  -l filename  | using external labels.\n";
  cerr << "  -w           | using weighted distance between faces.\n";
  cerr << "  -f           | do not use cache file.\n\n";

  cerr << "Unfolding\n";
  cerr << "  -ga filename | specify the ga config file. default = '"
      << default_config.ga_config_filename << "'\n";

  cerr << "  -lp filename | specify the LP config file. default = '"
      << default_config.lp_config_filename << "'\n";

  cerr << "  -s seed      | specify random seed\n";

  cerr << "  -r times     | retry times, default is "
      << default_config.max_retries << "\n";

  cerr << "  -weights fn  | using the specify weights to unfold the mesh.\n";

  cerr << "  -q           | quiet mode.\n";

  cerr << "  -bf          | specify the base face.\n";

  cerr << "  -rb          | random base face.\n";

  cerr << "  -pc          | using pixel checker for overlapping estimation.\n";

  cerr << "  -lc          | less cuts / don't cut flat edges.\n";

  cerr << "  -g           | disable GUI. dump outputs.\n";

  cerr << "  -tab         | add tabs in the net.\n";

  cerr << "  -nbb         | do not find best base face.\n\n";

  cerr <<"Net Optimization\n";
  cerr << "  -objective # | specify the objective for optimization\n";
  cerr << "     hull_area (minimize convex hull area)\n";
  cerr << "     cut_length (minimize cut edge length)\n\n";
  cerr << "     polygon (caging polygon, see -polygon)\n\n";

  cerr <<"Polygon Net Optimization\n";
  cerr << "  -polygon stencil_image.png/jpg | specify polygon for optimization\n";

  cerr << "Dumping SVG\n";
  cerr
      << "  -scale       | scale factor applied. both *.svg and *.ori will be affected.\n";
  cerr << "  -nl          | do not dump labels in SVG file\n";
  cerr << "  -lfs         | label font size scale [default=1.0]\n";
  cerr << "  -ecx         | extra cut x [default="
      << default_config.extra_cuts_x << "]\n";
  cerr << "  -ecy         | extra cut x [default="
      << default_config.extra_cuts_y << "]\n";

  cerr << "Others\n";
  cerr << "  -no_dump     | do not dump\n";

  cerr << endl;
  cerr << "-- Complied on " << __DATE__ << endl;
  cerr << "-- Report bugs to: Jyh-Ming Lien jmlien@cs.gmu.edu" << endl;
}

//-----------------------------------------------------------------------------

inline void run_ga_unfolder(Unfolder* unfolder)
{
  // using ga to find the best unfolding
  UnfoldingProblem problem(unfolder);

  if (problem.setup(config.ga_config_filename)) {
    // don't forget to set seed since it use another random generator instead of drand48 function
    problem.setSeed(config.seed);
    problem.run();
  } else {
    cerr << "Error! Failed to setup GA" << endl;
  }
}

inline void run_cluster_unfolder(Unfolder* unfolder)
{
  cout << "- Clustering module" << endl;
  // build K unfoldings and find associations between the faces on mesh
  // cluster faces that can folded together
  ClusterUnfolding problem(unfolder);

  if (problem.setup(config.cluster_config_filename)) {
    if (config.repair_mode && !config.label_filename.empty()) {
      problem.runRepair(config.label_filename);
    } else {
      // don't forget to set seed since it use another random generator instead of drand48 function
      //problem.setSeed(config.seed);
      problem.run();

      // Add all new unfolders here, no matter flattened or not...
      const auto new_unfodlers = problem.getUnfolders();
      cout << "- Done clustering, total unfolders = "
          << new_unfodlers.size() << endl;
      unfolders.insert(unfolders.end(), new_unfodlers.begin(),
          new_unfodlers.end());

      // To folded state
      updateUnfolding(0.0f);

      // Use random color for each cluster...
      randomColors();
    }
  } else {
    cerr << "Error! Failed to setup Clustering" << endl;
  }
} //run_cluster_unfolder()

bool readfromfiles()
{
  long vsize = 0;
  long fsize = 0;

  uint id = 0;

  auto start = clock();

  for (const auto& filename : filenames) {

    // use the same seed for all the component
    mathtool::srand48(config.seed);
    srand(config.seed);

    if(!config.quite)
      cout << "- [" << ++id << "/" << filenames.size() << "] Start reading "
           << filename << endl;

    config.filename = filename;

    model* m = new model();
    if (!m->build(config.filename, config.quite))
    {
      delete m;
      return false;
    }

    if(!config.quite)
      cout << "- Done reading " << m->v_size << " vertices and " << m->t_size
          << " facets" << endl;

    Unfolder* unfolder = new Unfolder(m, config);
    unfolder->measureModel();

    if (config.weights_filename.length() > 0) {
      // build from weights;
      cout<<"build from weights 1"<<endl;
      unfolder->buildFromWeights(config.weights_filename);
      cout<<"build from weights 2"<<endl;
      unfolder->rebuildModel();
      //optimize more from the current weights
      if (config.heuristic == CutHeuristic::GA){ run_ga_unfolder(unfolder); }
    } else if (!config.cuts_filename.empty()) {
      // build from cuts
      unfolder->buildFromCuts(config.cuts_filename);
      unfolder->rebuildModel();
    } else if (config.heuristic == CutHeuristic::CLUSTERING) {
      run_cluster_unfolder(unfolder);
    } else if (config.heuristic == CutHeuristic::LP) {
      //deprecated
      // LPUnfolding problem(unfolder);
      // showUnfold = false;
      //
      // if (problem.setup(config.lp_config_filename)) {
      //   // don't forget to set seed since it use another random generator instead of drand48 function
      //   //problem.setSeed(config.seed);
      //   problem.run();
      // } else {
      //   cerr << "Error! Failed to setup LP" << endl;
      // }
    } else if (config.heuristic == CutHeuristic::GA) {
      run_ga_unfolder(unfolder);
    } else if (!config.label_filename.empty()) {
      // exp mode:

      // show model
      showUnfold = false;

      // assign labels;
      auto labels = masc::util::readList<int>(config.label_filename);
      const auto m = unfolder->getModel();
      for (int i = 0; i < labels.size(); ++i) {
        m->tris[i].cluster_id = labels[i];
      }

    } else if (!config.score_filename.empty()) {
      // show model
      showUnfold = false;
      showScore = true;

      // assign labels;
      auto weights = masc::util::readList<float>(config.score_filename);
      const auto m = unfolder->getModel();
      for (int i = 0; i < weights.size(); ++i) {
        m->tris[i].score = weights[i];
      }

      vector<float> vs;
      // compute per vertex score
      for (int i = 0; i < m->v_size; ++i) {
        auto& v = m->vertices[i];
        auto total_area = 0.0;
        // weighted average
        for (auto fid : v.m_f) {
          v.score += m->tris[fid].area * m->tris[fid].score;
          total_area += m->tris[fid].area;
        }
        v.score /= total_area;
        vs.push_back(v.score);

        if (config.output_vertex_score) {
          masc::util::writeList(config.vertex_score_filename, vs);
        }
      }
    } else {
      // using basic heuristic
      unfolder->buildUnfolding();
    }

    cout << string(40, '-') << endl;

    if (unfolder->getCheckOverlappingCalls() > 0 && !config.quite) {
      cout << "- Total CO calls  = " << unfolder->getCheckOverlappingCalls()
           << endl;
      cout << "- Total CO time   = "
           << unfolder->getTotalCheckOverlappingTime() * 1.0 / CLOCKS_PER_SEC
           << " s" << endl;
      cout << "- Average CO time = "
           << unfolder->getTotalCheckOverlappingTime() * 1.0
              / unfolder->getCheckOverlappingCalls() / CLOCKS_PER_SEC << " s"
           << endl;
    }

    //cerr << "isFlattened = " << unfolder->isFlattened() << endl;

    if (unfolder->isFlattened()) {
      if (config.optimize_unfolding) {
        const vector<float> weights = unfolder->optimizeUnfolding();
        unfolder->buildFromWeights(weights);
      }

      if (config.find_best_base_face) {
        unfolder->findBestBaseFace();
      }
    }

    // In clustering mode, each cluster should be added to *unfolders*, skip here.
    if (config.heuristic != CutHeuristic::CLUSTERING)
    {
      unfolder->rebuildModel();
      unfolder->unfoldTo(0.0);
      unfolders.push_back(unfolder);
    }

    cout << string(40, '-') << endl;

  } //end for (const auto& filename)

  auto time_cost = (clock() - start) * 1.0 / CLOCKS_PER_SEC;

  int flattened = 0;

  for (auto unfolder : unfolders)
    if (unfolder->isFlattened())
      ++flattened;

  cout << "Total flattened = " << flattened << "/" << unfolders.size() << endl;
  cout << "Total time = " << time_cost << " s" << endl;

  computeCOM_R();

  return true;
}

void computeCOM_R() {
  //compute a bbox
  double box[6] = { FLT_MAX, -FLT_MAX, FLT_MAX, -FLT_MAX, FLT_MAX, -FLT_MAX };
  //-------------------------------------------------------------------------
  for (auto& u : unfolders) {
    auto& m = *u->getModel();
    for (int j = 0; j < m.v_size; j++) {
      Point3d& p = m.vertices[j].p;
      if (p[0] < box[0])
        box[0] = p[0];
      if (p[0] > box[1])
        box[1] = p[0];
      if (p[1] < box[2])
        box[2] = p[1];
      if (p[1] > box[3])
        box[3] = p[1];
      if (p[2] < box[4])
        box[4] = p[2];
      if (p[2] > box[5])
        box[5] = p[2];
    }				//j
  }				//i

  //-------------------------------------------------------------------------
  // compute center of mass and R...
  COM.set((box[1] + box[0]) / 2, (box[3] + box[2]) / 2, (box[5] + box[4]) / 2);

  //-------------------------------------------------------------------------
  R = 0;
  for (auto& u : unfolders) {
    auto& m = *u->getModel();
    for (int j = 0; j < m.v_size; j++) {
      Point3d& p = m.vertices[j].p;
      double d = (p - COM).normsqr();
      if (d > R)
        R = d;
    }    //j
  }    //i

  R = sqrt(R);
}

#endif //_BF_MINKOWSKI_SUM_H_
