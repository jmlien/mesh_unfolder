/*
 * Breeder.cpp
 *
 *  Created on: Mar 16, 2015
 *      Author: zhonghua
 */

#include "Breeder.h"
#include "Problem.h"

namespace masc {
namespace ga {

Breeder::Breeder(Problem* problem) {
  this->m_tournament_size = 2;
  this->m_problem = problem;
  this->m_children_size = 2;
}

Breeder::~Breeder() {
  // TODO Auto-generated destructor stub
}

bool Breeder::setup(const vector<string>& tokens) {
  for (int i = 1; i < tokens.size(); ++i) {
    const auto& token = tokens[i];

    if (token == "tournament-size") {
      this->m_tournament_size = stoi(tokens[++i]);
    } else if (token == "children-size") {
      this->m_children_size = stoi(tokens[++i]);
    } else {
      cerr << "- [GA] ! Error: Unknown token = " << token << endl;
      return false;
    }
  }
  return true;
}

void Breeder::print(ostream& out) const {
  out << "- [GA] breeder tournament-size=" << this->m_tournament_size
      << " children-size=" << this->m_children_size << endl;
}

vector<Individual*> Breeder::breed(const vector<Individual*>& population) {
  Random& random = this->m_problem->getRandom();

  vector<Individual*> children;

  while (children.size() < this->m_children_size) {
    // select two parents
    Individual* parent0 = this->select(population);
    Individual* parent1 = this->select(population);

    // copy from parents
    Individual* child0 = new Individual(*parent0);
    Individual* child1 = new Individual(*parent1);

    // do crossover
    child0->crossover(*child1, child0->getSpecies()->getCrossoverProb());

    // do mutation
    // mutation probability is tested on each gene
    child0->mutate();
    child1->mutate();

    children.push_back(child0);
    children.push_back(child1);
  }

  return children;
}

Individual* Breeder::select(const vector<Individual*>& population) {
  int init_index = drand48() * population.size();

  Individual* best = population[init_index];

  for (int i = 1; i < this->m_tournament_size; ++i) {
    int index = drand48() * population.size();
    Individual* ind = population[index];
    if (ind->getFitness() > best->getFitness())
      best = ind;
  }

  return best;
}

} /* namespace ga */
} /* namespace masc */
