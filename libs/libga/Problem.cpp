/*
 * Problem.cpp
 *
 *  Created on: Mar 16, 2015
 *      Author: zhonghua
 */

#include "Problem.h"

#include <fstream>
#include <cassert>
#include <ctime>
#include <cfloat>
#include <algorithm>
#include <iomanip>

#include "util/StreamHelper.h"

using namespace masc::ga::util;

namespace masc {
namespace ga {

Problem::Problem() {
  this->m_species = nullptr;
  this->m_breeder = nullptr;
  this->m_cur_generation = 0;
  this->m_max_generateions = 10;
  this->m_population_size = 10;
  this->m_goal_achieved = false;
  this->m_start_time = 0;
  this->m_generation_survived=0;

  this->m_best_ind.setFitness(-1e10);

}

Problem::~Problem() {
  delete this->m_species;
  this->m_species = nullptr;

  delete this->m_breeder;
  this->m_breeder = nullptr;

  for (auto ind : this->m_population) {
    delete ind;
  }

  this->m_population.clear();
}

bool Problem::setup(const string& filename) {
  ifstream in;
  in.open(filename, ios_base::in);
  if (!in.good()) {
    cerr << "- [GA] ! Error in Problem::setup. Cannot open file = " << filename
         << endl;
    return false;
  }

  // invoke implementation
  auto flag = this->setup(in);

  in.close();

  return flag;
}

bool Problem::setup(istream& in) {
  while (true) {
    vector<string> tokens = StreamHelper::tokenize(in);

    if (tokens.size() == 0)
      break;

    if (tokens[0] == "problem") {
      if (!this->setup(tokens)) {
        return false;
      }
    } else if (tokens[0] == "species") {
      // TODO use factory later
      this->m_species = new Species(this);

      if (!this->m_species->setup(tokens)) {
        return false;
      }
    } else if (tokens[0] == "breeder") {
      //TODO use factory later
      this->m_breeder = new Breeder(this);

      if (!this->m_breeder->setup(tokens)) {
        return false;
      }
    } else if (tokens[0] == "creators") {
      for (int i = 1; i < tokens.size(); ++i)
        this->m_creator_methods.push_back(tokens[i]);

    } else {
      return false;
    }
  }

  this->init();

  this->print(cout);

  // make sure species is created
  assert(m_species != nullptr);
  m_species->print(cout);

  // make sure breeder is created
  assert(m_breeder != nullptr);
  m_breeder->print(cout);

  return true;
}

bool Problem::setup(vector<string>& tokens) {
  this->setSeed(time(nullptr));

  for (int i = 0; i < tokens.size(); ++i) {
    const auto& token = tokens[i];

    if (token == "max-generation") {
      this->m_max_generateions = stoi(tokens[++i]);
    } else if (token == "pop-size") {
      this->m_population_size = stoi(tokens[++i]);
    } else if (token == "seed") {
      this->setSeed(stoi(tokens[++i]));
    }
  }
  return true;
}

// will be called once after setup,
// can be used to override parameters
void Problem::init() {
  //nothing to do here
}

void Problem::print(ostream& out) const
{
  out << "- [GA] problem max-generation=" << this->m_max_generateions << " pop-size="
      << this->m_population_size << endl;

  out << "- [GA] creators=";
  for (const string& creator : m_creator_methods)
    out << creator << " ";
  out << endl;
}

void Problem::run()
{
  m_start_time = clock();

  // generate initial population
  this->generatePopulation();

  for (auto i = 1; i <= this->m_max_generateions; ++i)
  {
    // breed children
    vector<Individual*> children = this->m_breeder->breed(this->m_population);

    // evaluate children
    for (Individual* child : children)
      this->evaluate(child);

    // update population
    this->mergePopulation(children);

    // callback
    this->generationDone(i);

    if(this->m_goal_achieved)
      break;
  }
  cerr << endl;
}

void Problem::generationDone(int generation)
{
  auto gen_best = this->getBestIndividual();

  if (gen_best->getFitness() > this->m_best_ind.getFitness()) {
    this->m_best_ind = *gen_best;
    this->m_generation_survived=0;
  }

  this->m_generation_survived++; //increase the number of survived generations

  auto avg_fitness = 0.0f;

  for (const auto& ind : m_population)
    avg_fitness += ind->getFitness();

  avg_fitness /= m_population.size();

  auto now = clock();

  auto time_in_sec = (now - m_start_time) * 1.0 / CLOCKS_PER_SEC;

  cerr << std::fixed << std::setprecision(4);

  cerr << "- [GA]  generation = " << generation << ", time = " << time_in_sec
       << " s, best = " << gen_best->getFitness() << "("<<m_generation_survived
       <<"), avg = " << avg_fitness<< ", "<<" m/c prob = "
       <<getSpecies()->getMutationProb()<<"/"<<getSpecies()->getCrossoverProb();
       //<<" cur="<<gen_best->getFitness();

  //this->m_gen_avg_fitness.push_back(avg_fitness);
  //this->m_gen_best_fitness.push_back(gen_best->getFitness());

  generationPrint(cerr, *gen_best);

  cerr << "\r" << flush;
  //cout << endl;
}

void Problem::generationPrint(ostream& out, const Individual& gen_best) {
  // do nothing here
}

Individual* Problem::generateIndividual() {
  return this->m_species->createNewIndvidual();
}

void Problem::mergePopulation(const vector<Individual*>& children)
{
  // sort the population based on fitness
  sort(this->m_population.begin(), this->m_population.end(),
      [](const Individual* a, const Individual* b) {
        return *a < *b;
      });

  for (auto i = 0u; i < children.size(); ++i) {
    // kill parent
    delete m_population[i];
    // put child in population
    m_population[i] = children[i];
  }
}

void Problem::generatePopulation()
{
  cout << "- [GA]  generating population..." << endl;

  for (int i = 0; i < this->m_population_size; ++i) {
    auto ind = this->generateIndividual();

    this->evaluate(ind);

    this->m_population.push_back(ind);

    if (ind->getFitness() > this->m_best_ind.getFitness()) {
      this->m_best_ind = *ind;
      cout << "- [GA]  [" << i << "/" << m_population_size << "] best so far = "
           << m_best_ind.getFitness() << "\r"<<flush;
    }
  }

  cout << "\n- [GA]  finished" << endl;
}

Individual* Problem::getBestIndividual()
{
  Individual* best = m_population[0];

  for (int i = 0; i < m_population.size(); ++i) {
    Individual* ind = m_population[i];
    if (ind->getFitness() > best->getFitness())
      best = ind;
  }

  return best;
}

Individual* Problem::getWorstIndividual(int& index) {
  Individual* worst = m_population[0];
  index = 0;

  for (int i = 1; i < m_population.size(); ++i) {
    Individual* ind = m_population[i];
    if (ind->getFitness() < worst->getFitness()) {
      worst = ind;
      index = i;
    }
  }

  return worst;
}

Species* Problem::getSpecies() {
  return this->m_species;
}

//////////////////////////////////////////////////
// access
//////////////////////////////////////////////////

void Problem::setSeed(int seed) {
  this->m_random.setSeed(seed);
}

Random& Problem::getRandom() {
  return this->m_random;
}

} /* namespace ga */
} /* namespace masc */
