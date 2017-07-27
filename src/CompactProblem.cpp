/*
 * CompactProblem.cpp
 *
 *  Created on: Apr 13, 2016
 *      Author: zxi
 */

#include "CompactProblem.h"

#include <cfloat>

#include "unfolder.h"
#include "mathtool/Box.h"

namespace masc {
namespace unfolding {

//
class DiscreteAssigner: public FoldingAngleAssigner {
public:
  DiscreteAssigner(Unfolder* unfolder) :
      FoldingAngleAssigner(unfolder) {
  }

  virtual ~DiscreteAssigner() {
  }

  // assign folding angles for an individual
  virtual void assign(Individual* ind) override {
    for (auto i = 0; i < ind->getGenome().size(); ++i) {
      float gene = 0.0f;
      switch (mathtool::lrand48() % 4) {
      case 0:
        gene = 0.5f; // folding angle = 0
        break;
      case 1:
        gene = 0.0f; // folding angle = -PI
        break;
      case 2:
        gene = 1.0f;  // folding angle = PI
        break;
      case 3:
        gene = this->m_unfolder->getInitFoldingAngles()[i] / PI / 2 + 0.5;
        break;
      }
      ind->setGene(i, gene);
    }
  }
};

class FlatAssigner: public FoldingAngleAssigner {
public:
  FlatAssigner(Unfolder* unfolder) :
      FoldingAngleAssigner(unfolder) {
  }

  virtual ~FlatAssigner() {
  }

  // assign folding angles for an individual
  virtual void assign(Individual* ind) override {
    for (auto i = 0; i < ind->getGenome().size(); ++i) {
      float gene = 0.0f;
      switch (mathtool::lrand48() % 3) {
      case 0:
        gene = 0.5f; // folding angle = 0
        break;
      case 1:
        gene = 0.0f; // folding angle = -PI
        break;
      case 2:
        gene = 1.0f;  // folding angle = PI
        break;
      }
      ind->setGene(i, gene);
    }
  }
};

// folding angle is from [0, target_folding_angle] or [target_folding_angle, 0]
class InrangeAssigner: public FoldingAngleAssigner {
public:
  InrangeAssigner(Unfolder* unfolder) :
      FoldingAngleAssigner(unfolder) {
  }

  virtual ~InrangeAssigner() {

  }

  // assign folding angles for an individual
  virtual void assign(Individual* ind) override {
    for (auto i = 0; i < ind->getGenome().size(); ++i) {
      float gene = this->m_unfolder->getInitFoldingAngles()[i]
          * mathtool::drand48() / PI / 2 + 0.5;
      ind->setGene(i, gene);
    }
  }
};

CompactProblem::CompactProblem(Unfolder* unfolder) {

  this->m_unfolder = unfolder;
  this->m_best_compactness = 0;
  this->m_last_compactness = 0;

  this->m_assigners.push_back(
      std::move(
          std::unique_ptr<FoldingAngleAssigner>(
              new DiscreteAssigner(unfolder))));

  this->m_assigners.push_back(
      std::move(
          std::unique_ptr<FoldingAngleAssigner>(
              new InrangeAssigner(unfolder))));

  this->m_assigners.push_back(
      std::move(
          std::unique_ptr<FoldingAngleAssigner>(
              new FlatAssigner(unfolder))));
}

CompactProblem::~CompactProblem() {
  this->m_unfolder = nullptr;
}

///////////////////////////////////////////////////////////////////////////////
// Protected
///////////////////////////////////////////////////////////////////////////////

void CompactProblem::init() {
  // only mutate fold edges...
  this->m_species->setGenomeSize(this->m_unfolder->getModel()->t_size - 1);
}

Individual* CompactProblem::generateIndividual() {
  // by default use random genome
  auto ind = Problem::generateIndividual();

  int index = this->getRandom().nextDoubleUniform() * this->m_assigners.size();

  if (index <= this->m_assigners.size() - 1) {
    this->m_assigners.at(index)->assign(ind);
  }

  return ind;
}

vector<double> CompactProblem::genome2foldingAngles(
    const vector<float>& genome) {
  vector<double> fold_angles(genome.size());

  // map 0 ~ 1 -> -PI ~ PI
  for (int i = 0; i < genome.size(); ++i) {
    fold_angles[i] = (genome[i] - 0.5) * 2 * PI;
    if (genome[i] < 1e-6)
      fold_angles[i] = -PI;
    if (genome[i] > (1.0 - 1e-6))
      fold_angles[i] = PI;
  }

  auto full = this->m_unfolder->getFullFoldingAngles(fold_angles);

  return full;
}

void CompactProblem::evaluate(Individual* ind) {
  const auto& config = this->m_unfolder->getConfig();

  auto full = this->genome2foldingAngles(ind->getGenome());

  this->m_unfolder->unfoldTo(full);

  auto fitness = this->evaulte();

  ind->setFitness(fitness);

  if (this->m_best_ind.getFitness() > m_best_compactness * 1.001) {
    m_best_compactness = this->m_best_ind.getFitness();
    auto full = this->genome2foldingAngles(this->m_best_ind.getGenome());
    m_evolving_seqs.push_back(full);
  }
}

float CompactProblem::evaulte() {
  if (this->m_unfolder->checkCollision())
    return 0.0f;

  mathtool::Box3d box(this->m_unfolder->getVertices());

  return 1.0f / box.getSumOfDim();
}

void CompactProblem::generationDone(int generation) {
  const auto config = this->m_unfolder->getConfig();

  Problem::generationDone(generation);

  // finished
  if (this->m_best_ind.getFitness() == 1.0
      || generation == this->m_max_generateions) {

    this->m_unfolder->unfoldTo(0.0);
    cout << " - original mesh compactness = " << this->evaulte() << endl;

    this->m_unfolder->unfoldTo(1.0);
    cout << " - net compactness = " << this->evaulte() << endl;

    // fold to best
    auto full = this->genome2foldingAngles(m_best_ind.getGenome());

    this->m_unfolder->unfoldTo(full);
    cout << " - best compactness = " << this->evaulte() << endl;
  }

}
} /* namespace unfolding */
} /* namespace masc */
