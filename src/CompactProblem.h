/*
 * CompactProblem.h
 *
 *  Created on: Apr 13, 2016
 *      Author: zxi
 */

#ifndef SRC_COMPACTPROBLEM_H_
#define SRC_COMPACTPROBLEM_H_

#include <memory>
#include "libga/Problem.h"
using masc::ga::Individual;

class Unfolder;
class FoldingAngleAssigner;

namespace masc {
namespace unfolding {

class FoldingAngleAssigner {
public:
  // assign folding angles for an individual
  virtual void assign(Individual* ind) = 0;
  virtual ~FoldingAngleAssigner() {
    this->m_unfolder = nullptr;
  }
protected:
  FoldingAngleAssigner(Unfolder* unfolder) {
    this->m_unfolder = unfolder;
  }

  Unfolder* m_unfolder;
};

class CompactProblem: public ga::Problem {
public:
  CompactProblem(Unfolder* unfolder);
  virtual ~CompactProblem();

  virtual void evaluate(Individual* ind) override;

  virtual Individual* generateIndividual() override;

  const vector<vector<double>>& getEvolvingSeqs() {
    return m_evolving_seqs;
  }

protected:
  virtual void init() override;
  virtual void generationDone(int generation) override;

  vector<double> genome2foldingAngles(const vector<float>& genome);

  /// evaluated the compactness of the folded net
  float evaulte();

  // does not own object
  Unfolder* m_unfolder;

  // best compactness
  double m_best_compactness;

  // last best one
  double m_last_compactness;

  // sequence of "best" folding angles
  vector<vector<double>> m_evolving_seqs;

  vector<unique_ptr<FoldingAngleAssigner>> m_assigners;
};

} /* namespace unfolding */
} /* namespace masc */

#endif /* SRC_COMPACTPROBLEM_H_ */
