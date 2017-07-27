/*
 * Breeder.h
 *
 *  Created on: Mar 16, 2015
 *      Author: zhonghua
 */

#ifndef BREEDER_H_
#define BREEDER_H_

#include "Species.h"

namespace masc {
namespace ga {


// forward declaration
class Problem;

/*
 * currently, breeder is simply implemented using tournamentSelection to select
 * two individuals and then crossover and mutate which can be improved/extended later
 */
class Breeder {
public:
    Breeder(Problem* problem);
    virtual ~Breeder();
    virtual bool setup(const vector<string>& tokens);

    // breed new individuals
    virtual vector<Individual*> breed(const vector<Individual*>& population);

    // select the best one from population
    virtual Individual* select(const vector<Individual*>& population);

    virtual void print(ostream& out) const;

protected:

    // how many children to breed each generation
    int m_children_size;

    // among how many parents select the best one
    int m_tournament_size;

    // pointer to the problem
    Problem* m_problem;
};

} /* namespace ga */
} /* namespace masc */

#endif /* BREEDER_H_ */
