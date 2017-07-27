/*
 * Individual.h
 *
 *  Created on: Mar 16, 2015
 *      Author: zhonghua
 */

#ifndef INDIVIDUAL_H_
#define INDIVIDUAL_H_

#include <vector>
using namespace std;

namespace masc {
namespace ga {

class Species;

class Individual {
public:
  Individual();
  Individual(int genome_size, Species* species);
  Individual(const Individual& ind);
  virtual ~Individual();

  ////////////////////////////////////////////////////

  // reset the individual
  virtual void reset();

  // mutate the individual
  virtual void mutate();

  // crossover two individuals
  virtual void crossover(Individual& ind, float crossover_prob);

  // compare two individual based on their fitness
  bool operator<(const Individual& other) const;

  ////////////////////////////////////////////////////
  // access
  ////////////////////////////////////////////////////

  void setGenome(const vector<float>& genome);
  const vector<float>& getGenome() const;

  void setGene(const int index, const float gene);
  const float getGene(const int index) const;

  const bool getEvaluated() const;
  void setEvaluated(const bool evaluated);

  const double getFitness() const;
  void setFitness(const double fitness);

  const bool getIsValid() const;
  void setIsValid(const bool is_valid);

  Species* getSpecies();

  // Species can access genome directly
  friend class Species;

protected:
  const float getRandomGene() const;

private:

  Species* m_species;

  vector<float> m_genome;

  // number of genes
  int m_genome_size;

  // whether this individual has been evaluated or not
  bool m_evaluated;

  // the fitness score of the individual
  // smaller means the worse
  // bigger means the better
  double m_fitness;

  /// whether the individual is valid or not.
  bool m_is_valid;
};

} /* namespace ga */
} /* namespace masc */

#endif /* INDIVIDUAL_H_ */
