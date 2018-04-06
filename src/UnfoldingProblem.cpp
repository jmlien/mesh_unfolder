/*
 * UnfoldingProblem.cpp
 *
 *  Created on: Mar 17, 2015
 *      Author: zhonghua
 */

#include "UnfoldingProblem.h"
#include "unfolder.h"
#include "UnfoldingEvaluator.h"
#include "Splitter.h"
#include "OverlappingChecker.h"

namespace masc {
namespace unfolding {

///////////////////////////////////////////////////////////////////////////////
// Unfolding Problem
///////////////////////////////////////////////////////////////////////////////
UnfoldingProblem::UnfoldingProblem(Unfolder* unfolder) {

  this->m_unfolder = unfolder;
  //this->m_checker.reset(new PixelChecker(unfolder));
  //this->m_checker.reset(new ITreeChecker(unfolder));
  this->m_checker.reset(new BruteForceChecker(unfolder));

  if (unfolder->getConfig().pixel_checker) {
    this->m_evaluator.reset(new AreaEvaluator(unfolder));
  } else if (loadLearningParams(unfolder->getConfig().params_file_path)) {
    this->m_evaluator.reset(new LearningEvaluator(unfolder, &m_params));
  } else {
    this->m_evaluator.reset(new OverlappingEvaluator(unfolder));
  }

  switch (unfolder->getConfig().opt_obj) {
  case Objective::CUT_LENGTH:
    this->m_net_evaluator.reset(new CutLengthEvaluator);
    break;
  case Objective::HULL_AREA:
    this->m_net_evaluator.reset(new HullAreaEvaluator);
    break;
  case Objective::POLYGON_FIT:
    this->m_net_evaluator.reset(new PolygonFitEvaluator(unfolder->getConfig().stencil_filename));
    break;
  default:
    break;
  }

  this->m_best_ratio = 1e10;
}

UnfoldingProblem::~UnfoldingProblem() {
  this->m_unfolder = nullptr;
  m_spliiters.clear();
}

void UnfoldingProblem::init() {
  // genome size = edge size
  this->m_species->setGenomeSize(this->m_unfolder->getModel()->e_size);

  for (const string& creator : this->m_creator_methods) {
    this->m_spliiters.push_back(
        std::move(
            std::unique_ptr<Splitter>(Splitter::createSplitter(creator))));
  }

  for (const auto& splitter : m_spliiters)
    splitter->measure(this->m_unfolder->getModel());

  return;

  // TODO(zxi) one time creators...
#if 0
  this->m_oneshot_spliiters.push_back(
      std::move(std::unique_ptr<Splitter>(new MinimumPerimeterSplitter())));
  this->m_oneshot_spliiters.push_back(
      std::move(std::unique_ptr<Splitter>(new MaximumPerimeterSplitter())));

  for (const auto& splitter : m_oneshot_spliiters) {
    splitter->measure(this->m_unfolder->getModel());

    auto indPtr = Problem::generateIndividual();
    auto weights = splitter->assignWeights(this->m_unfolder->getModel(),
        this->m_unfolder->getConfig());

    indPtr->setGenome(weights);
    this->evaluate(indPtr);
    cerr << "fitness = " << indPtr->getFitness() << endl;
    this->m_population.push_back(indPtr);
  }
#endif

}

void UnfoldingProblem::run() {

  auto start_time = clock();
  auto& config = this->m_unfolder->getConfig();
  this->m_max_generateions=min(this->m_max_generateions, config.ga_max_gen);

  //there is already a weight, so this is a rerun...
  if(this->m_unfolder->getEdgeWeights().empty()==false) //this model is build from existing weights
  {
    cout<<" - added weight as the best individual before creating more..."<<endl;
    auto ind = Problem::generateIndividual();
    ind->setGenome(this->m_unfolder->getEdgeWeights());
    this->evaluate(ind);
    this->m_population.push_back(ind);
    this->m_best_ind = *ind;
  }

  // generate initial population
  this->generatePopulation();

  int generation=1;
  while(m_goal_achieved==false)
  {
    for (; generation <= this->m_max_generateions; ++generation)
    {
      // breed children
      vector<Individual*> children = this->m_breeder->breed(this->m_population);

      // evaluate children
      for (Individual* child : children)
        this->evaluate(child);

      // update population
      this->mergePopulation(children);

      // callback
      this->generationDone(generation);

      if (m_goal_achieved)
        break;
    }//end for

    if(this->m_best_ind.getIsValid()==false)//failed
    {
      m_goal_achieved=false;
      cerr<<"! Unfolding failed. Would you like to try more interations (say \"no\" or provide a number)? ";
      string answer;
      cin>>answer;
      if(answer=="no" || answer=="No" || answer=="NO"){
        break;
      }
      this->m_max_generateions+=atoi(answer.c_str());
    }//end m_goal_achieved

  }//end while

  cerr << endl;
}

Individual* UnfoldingProblem::generateIndividual() {
  auto indPtr = Problem::generateIndividual();

  // use a random splitter to generate weights
  const auto& splitter = this->m_spliiters.at(
      this->m_spliiters.size() * this->getRandom().nextDoubleUniform());

  auto weights = splitter->assignWeights(this->m_unfolder->getModel(),
      this->m_unfolder->getConfig());

  indPtr->setGenome(weights);

  return indPtr;
}

void UnfoldingProblem::trimGene(Individual* ind) {
  auto m = this->m_unfolder->getModel();
  for (int i = 0; i < m->e_size; ++i) {
    const edge& e = m->edges[i];

    // won't accept negative weight, will all become zero.
    if (e.type == 'd') {
      ind->setGene(i, EdgeWeight::DiagnalEdge);
    } else if (e.type == 'p') {
      ind->setGene(i, EdgeWeight::FlatEdge);
    }
  }
}

void UnfoldingProblem::evaluate(Individual* ind) {
  const auto& config = this->m_unfolder->getConfig();
// avoid cutting on flat edges e.g. diagonals
  if (config.less_cuts) {
    this->trimGene(ind);
  }

  const auto fitness = this->m_evaluator->evaluate(ind);
  ind->setFitness(fitness);

  if (fitness < 0 && m_unfolder->getLastOverlapCount() < config.ga_local_adj) {
    //TODO(zxi) do some local adjustment
    // collect edges of overlapping triangles.
    set<uint> overlap_edges;
    for (const triangle& f : this->m_unfolder->getModel()->tris) {
      if (!f.overlapped)
        continue;
      for (int i = 0; i < 3; ++i)
        overlap_edges.insert(f.e[i]);
    }

    for (const uint eid : overlap_edges) {
      if (this->getRandom().nextDoubleUniform() > 0.5)
        continue;

      // Cut the edge
      ind->setGene(eid, 1.0);
    }

    // Re-evaluate.
    const auto fitness = this->m_evaluator->evaluate(ind);
    ind->setFitness(fitness);
  }

  // Post evaluate the net
  if (ind->getFitness() == 0 && this->m_net_evaluator != nullptr) {
    double post_fitness = this->m_net_evaluator->evaluate(this->m_unfolder);
    assert(post_fitness>=0);
    ind->setFitness(post_fitness);
  }

// Evolved a better one, record it!
  if (this->m_best_ind < *ind) {
    UnfoldingState state(ind->getGenome());
    this->m_evolving_seqs.push_back(state);
  }

// is a valid solution
  if (this->m_unfolder->isFlattened()) {
    ind->setIsValid(true);
  }
}

void UnfoldingProblem::generationDone(int generation) {
  const auto& config = this->m_unfolder->getConfig();

  Problem::generationDone(generation);

// finished
  if ((this->m_best_ind.getIsValid() && config.early_stop)
      || generation == this->m_max_generateions) {

    // Update the evolving sequences
    std::cerr << endl;
    std::cerr << "- Total evolving sequences: " << this->m_evolving_seqs.size()
        << std::endl;
    this->m_unfolder->setEvolvingSeqs(this->m_evolving_seqs);

    // rebuild unfolding from the best individual
    this->m_unfolder->buildFromWeights(m_best_ind.getGenome());
    this->m_unfolder->rebuildModel();

    // goal achieved
    this->m_goal_achieved = true;

//    for (int i = 0; i < m_gen_avg_fitness.size(); ++i)
//      cout << "gen " << i << " avg_fitness = " << m_gen_avg_fitness[i]
//          << " best_finess = " << m_gen_best_fitness[i] << endl;
  }
}

void UnfoldingProblem::generationPrint(ostream& out,
    const Individual& gen_best) {
  if (gen_best.getIsValid()) {
    this->m_unfolder->buildFromWeights(gen_best.getGenome(), false);
    float obj_value = 1.0 / gen_best.getFitness();
    out << " objective = " << obj_value;
  }
}

bool UnfoldingProblem::loadLearningParams(const string& inputsFile) {
  if (inputsFile == "")
    return false;

  istream *input;
  ifstream infile;
  istringstream inString;

  infile.open(inputsFile.c_str(), ifstream::in);

  if (!infile) {
    cerr << "Could not open learning parameter file " << inputsFile << endl;
    return false;
  }

  input = &(infile);

  string name;
  bool fBlockComment = false;
  while (!input->eof()) {

    // Skip comments and empty lines
    std::string str;
    std::getline(*input, str);
    if (str.length() >= 2 && str.substr(0, 2) == "/*") {
      fBlockComment = true;
    } else if (str == "*/") {
      fBlockComment = false;
    }
    if (fBlockComment || str == "" || str[0] == '#') {
      continue;
    }

    // otherwise parse strings
    stringstream s(str);
    std::string key;
    std::string value;
    std::getline(s, key, '\t');      //read thru tab
    std::getline(s, value);          //read thru newline
    if (value.empty()) {
      continue;
    }
    m_params[key] = value;
    std::cout << "Loaded learning paramter " << key << " = "
        << atof(value.c_str()) << std::endl;
  }

  infile.close();
  return true;
}

} /* namespace unfolding */
} /* namespace masc */
