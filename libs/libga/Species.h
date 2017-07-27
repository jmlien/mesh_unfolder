/*
 * Species.h
 *
 *  Created on: Mar 16, 2015
 *      Author: zhonghua
 */

#ifndef SPECIES_H_
#define SPECIES_H_

#include <iostream>
#include <string>
#include <vector>
#include <cassert>
using namespace std;

#include "Individual.h"

namespace masc {
namespace ga {

class Problem;

enum class CrossoverType {
  ONE_POINT, TWO_POINTS, UNIFORM
};

enum class MutationType {
  RESET, GAUSSIAN
};

class CrossoverTypeConverter {
public:
  static CrossoverType parse(const string& type) {
    if (type == toString(CrossoverType::ONE_POINT)) {
      return CrossoverType::ONE_POINT;
    } else if (type == toString(CrossoverType::TWO_POINTS)) {
      return CrossoverType::TWO_POINTS;
    } else if (type == toString(CrossoverType::UNIFORM)) {
      return CrossoverType::UNIFORM;
    } else {
      cerr << "!Error. Unknown crossover type = " << type << endl;
      assert(false);
    }

    return CrossoverType::ONE_POINT;
  }

  static string toString(const CrossoverType& type) {
    switch (type) {
    case CrossoverType::ONE_POINT: {
      return "one-point";
      break;
    }
    case CrossoverType::TWO_POINTS: {
      return "two-points";
      break;
    }
    case CrossoverType::UNIFORM: {
      return "uniform";
      break;
    }
    default: {
      cerr << "!Error. Unknown crossover type = " << (int) type << endl;
      break;
    }
    }

    return "unknown";
  }
};

class MutationTypeConverter {
public:
  static const MutationType parse(const string& type) {
    if (type == toString(MutationType::RESET)) {
      return MutationType::RESET;
    } else if (type == toString(MutationType::GAUSSIAN)) {
      return MutationType::GAUSSIAN;
    } else {
      cerr << "!Error. Unknown crossover type = " << type << endl;
      assert(false);
    }

    return MutationType::RESET;

  }

  static const string toString(const MutationType& type) {
    switch (type) {
    case MutationType::RESET: {
      return "reset";
      break;
    }
    case MutationType::GAUSSIAN: {
      return "gaussian";
      break;
    }
    default: {
      cerr << "!Error. Unknown mutation type = " << (int) type << endl;
      break;
    }
    }

    return "unknown";
  }
};

class Species {
public:
  Species(Problem* problem);
  virtual ~Species();

  // setup from tokens
  virtual bool setup(const vector<string>& tokens);

  // create a new individual to fill the population
  virtual Individual* createNewIndvidual();

  // print species info
  virtual void print(ostream& out);

  ///////////////////////////////////
  // access
  ///////////////////////////////////
  void setMinGene(const float min_gene);
  const float getMinGene() const;

  void setMaxGene(const float max_gene);
  const float getMaxGene() const;

  void setGenomeSize(int genome_size);
  const int getGenomeSize() const;

  void setChunkSize(int chunk_size);
  const int getChunkSize() const;

  void setMutationProb(const float prob);
  const float getMutationProb() const;

  void setMutationGaussianStddev(const double stddev);
  const double getMutationGaussianStddev() const;

  void setCrossoverProb(const float prob);
  const float getCrossoverProb() const;

  const MutationType getMutationType() const;
  const CrossoverType getCrossoverType() const;

  Problem* getProblem();

protected:

  float m_min_gene;
  float m_max_gene;

  int m_genome_size;
  int m_chunk_size;

  float m_mutation_prob;
  float m_crossover_prob;

  CrossoverType m_crossover_type;
  MutationType m_mutation_type;

  // standard deviation of gaussian mutation
  double m_mutation_gaussian_stddev;

  Problem* m_problem;

};

} /* namespace ga */
} /* namespace masc */

#endif /* SPECIES_H_ */
