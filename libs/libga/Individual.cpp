/*
 * Individual.cpp
 *
 *  Created on: Mar 16, 2015
 *      Author: zhonghua
 */

#include "Individual.h"

#include <cfloat>
#include "Species.h"
#include "Problem.h"

namespace masc {
namespace ga {

Individual::Individual() {

  this->m_evaluated = false;
  this->m_fitness = DBL_MIN;
  this->m_genome_size = 1;

  this->m_species = nullptr;

  this->m_is_valid = false;
}

Individual::Individual(int genome_size, Species* species) {

  this->m_evaluated = false;
  this->m_fitness = DBL_MIN;
  this->m_genome_size = genome_size;

  this->m_genome = vector<float>(genome_size);
  this->m_species = species;

  this->m_is_valid = false;

  this->reset();
}

// create new copy of individual
Individual::Individual(const Individual& ind) {
  this->m_evaluated = false;
  this->m_fitness = DBL_MIN;
  this->m_is_valid = false;
  this->m_genome_size = ind.m_genome_size;
  this->m_genome = ind.m_genome;
  this->m_species = ind.m_species;
}

Individual::~Individual() {
  this->m_species = nullptr;
}

void Individual::reset() {
  for (auto& gene : this->m_genome) {
    gene = this->getRandomGene();
  }

}

void Individual::mutate() {
  Random& random = this->m_species->getProblem()->getRandom();

  for (int i = 0; i < this->m_genome_size; ++i) {
    if (random.nextDoubleUniform() > this->m_species->getMutationProb())
      continue;

    // mutation happens

    float mutated_gene = this->getGene(i);

    switch (this->m_species->getMutationType()) {
    case MutationType::RESET: {
      mutated_gene = this->getRandomGene();
      break;
    }
    case MutationType::GAUSSIAN: {
      do {
        mutated_gene = random.nextDoubleGuasian(this->getGene(i),
            this->m_species->getMutationGaussianStddev());
      } while (mutated_gene > this->m_species->getMaxGene()
          || mutated_gene < this->m_species->getMinGene());
      break;
    }
    default: {
      assert(false);
      break;
    }
    }

    this->setGene(i, mutated_gene);
  }

  this->m_evaluated = false;
}

void Individual::crossover(Individual& ind, float crossover_prob) {
  Random& random = this->m_species->getProblem()->getRandom();
  const auto chunk_size = this->m_species->getChunkSize();

  switch (this->m_species->getCrossoverType()) {
  case CrossoverType::ONE_POINT: {
    // [0, genome_size-1]
    int point = random.nextDoubleUniform() * (this->m_genome_size / chunk_size);

    // swap [0~point] chunks
    for (int x = 0; x < point * chunk_size; ++x) {
      if (random.nextDoubleUniform() < crossover_prob)
        swap(this->m_genome[x], ind.m_genome[x]);
    }

    break;
  }
  case CrossoverType::TWO_POINTS: {
    int point1 = random.nextDoubleUniform()
        * (this->m_genome_size / chunk_size);
    int point2 = random.nextDoubleUniform()
        * (this->m_genome_size / chunk_size);

    if (point1 > point2)
      swap(point1, point2);

    // swap [point1~point2] chunks
    for (int x = point1 * chunk_size; x < point2 * chunk_size; ++x) {
      if (random.nextDoubleUniform() < crossover_prob)
        swap(this->m_genome[x], ind.m_genome[x]);
    }

    break;
  }
  case CrossoverType::UNIFORM: {
    for (int x = 0; x < this->m_genome_size; ++x) {
      if (random.nextDoubleUniform() < crossover_prob)
        swap(this->m_genome[x], ind.m_genome[x]);
    }
    break;
  }
  default: {
    assert(false);
    break;
  }
  }
  // new individuals need to be evaluated
  this->m_evaluated = false;
  ind.m_evaluated = false;
}

//////////////////////////////////////////////////

const float Individual::getRandomGene() const {
  Random& random = this->m_species->getProblem()->getRandom();

  return random.nextDoubleUniform()
      * (this->m_species->getMaxGene() - this->m_species->getMinGene())
      + this->m_species->getMinGene();
}

bool Individual::operator<(const Individual& other) const {
  return this->m_fitness < other.getFitness();
}

///////////////////////////////////////////////
// access
///////////////////////////////////////////////

void Individual::setGenome(const vector<float>& genome) {
  this->m_genome = genome;
}

const vector<float>& Individual::getGenome() const {
  return this->m_genome;
}

const bool Individual::getEvaluated() const {
  return this->m_evaluated;
}

void Individual::setEvaluated(const bool evaluated) {
  this->m_evaluated = evaluated;
}

const double Individual::getFitness() const {
  return this->m_fitness;
}

void Individual::setFitness(const double fitness) {
  // set the evaluated property
  this->setEvaluated(true);

  this->m_fitness = fitness;
}

const bool Individual::getIsValid() const {
  return this->m_is_valid;
}

void Individual::setIsValid(const bool is_valid) {
  this->m_is_valid = is_valid;
}

void Individual::setGene(const int index, const float gene) {
  if (gene > this->m_species->getMaxGene())
    this->m_genome.at(index) = this->m_species->getMaxGene();
  else if (gene < this->m_species->getMinGene())
    this->m_genome.at(index) = this->m_species->getMinGene();
  else
    this->m_genome.at(index) = gene;
}

const float Individual::getGene(const int index) const {
  return this->m_genome.at(index);
}

Species* Individual::getSpecies() {
  return this->m_species;
}

} /* namespace ga */
} /* namespace masc */
