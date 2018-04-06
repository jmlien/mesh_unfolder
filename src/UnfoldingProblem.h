/*
 * UnfoldingProblem.h
 *
 *  Created on: Mar 17, 2015
 *      Author: zhonghua
 */

#ifndef UNFOLDINGPROBLEM_H_
#define UNFOLDINGPROBLEM_H_

#include "libga/Problem.h"
#include "UnfoldingState.h"
#include <memory>
#include <map>

class Unfolder;

using namespace masc::ga;

namespace masc {
class Splitter;
namespace unfolding {

class OverlappingChecker;
class UnfoldingEvaluator;
class NetEvaluator;

///////////////////////////////////////////////////////////////////////////////
// UnfoldingProblem
///////////////////////////////////////////////////////////////////////////////

class UnfoldingProblem: public Problem {
public:
  UnfoldingProblem(Unfolder* unfolder);
  virtual ~UnfoldingProblem();

  // evaluate an individual
  void evaluate(Individual* ind) override;

  // Generate an individual
  Individual* generateIndividual() override;

  // callback function after each generation for print
  void generationPrint(ostream& out, const Individual& gen_best) override;

  // override to ask user input 
  void run() override;

protected:
  // override to set genome size on the fly
  void init() override;

  void generationDone(int generation) override;

  bool loadLearningParams(const string& inputFile);

  // trim gene
  // 1. avoid to cut diagonal edges
  void trimGene(Individual* ind);

  // does not own object
  Unfolder* m_unfolder;

  double m_best_ratio;

  std::unique_ptr<UnfoldingEvaluator> m_evaluator;
  std::unique_ptr<NetEvaluator> m_net_evaluator;

  std::unique_ptr<OverlappingChecker> m_checker;

  std::vector<std::unique_ptr<Splitter>> m_spliiters;
  std::vector<std::unique_ptr<Splitter>> m_oneshot_spliiters;

  // best states during evolving
  std::vector<UnfoldingState> m_evolving_seqs;

  // parameters for learning
  std::map<string,string> m_params;
};

} /* namespace unfolding */
} /* namespace masc */

#endif /* UNFOLDINGPROBLEM_H_ */
