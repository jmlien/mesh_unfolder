/*
 * Problem.h
 *
 *  Created on: Mar 16, 2015
 *      Author: zhonghua
 */

#ifndef PROBLEM_H_
#define PROBLEM_H_

#include <string>
#include <iostream>
#include <vector>
#include <ctime>
#include <random>
using namespace std;

#include "Species.h"
#include "Individual.h"
#include "Breeder.h"

#ifdef _WIN32
#include "rand48.h"
#endif

#include "util/Random.h"
using namespace masc::ga::util;

namespace masc {
namespace ga {

class Problem {
public:
    Problem();
    virtual ~Problem();

    bool setup(const string& filename);

    //////////////////////////////////////////
    // ga
    //////////////////////////////////////////

    virtual void run();

    virtual Individual* generateIndividual();

    // evaluate an individual
    virtual void evaluate(Individual* ind) = 0;

    // generate random individuals
    void generatePopulation();

    // merge children with population of last generation
    void mergePopulation(const vector<Individual*>& children);

    //////////////////////////////////////////

    virtual void print(ostream& out) const;

    //////////////////////////////////////////
    // access
    //////////////////////////////////////////

    void setSeed(int seed);
    Random& getRandom();


    Individual* getBestIndividual();
    Individual* getWorstIndividual(int& index);

    Species* getSpecies();

    const int getPopulationSize() const;
    const int getCurGeneration() const;
    const int getMaxGenerations() const;

	  bool isGoalAchieved() const { return m_goal_achieved; }

protected:

    // will be called once after setup,
    // can be used to override parameters
    void virtual init();

    // setup the entire problem from stream
    virtual bool setup(istream& in);

    // setup the problem from tokens
    virtual bool setup(vector<string>& tokens);

    // callback function after each generation
    virtual void generationDone(int generation);

    // callback function after each generation for print
    virtual void generationPrint(ostream& out, const Individual& gen_best);


    // species of the problem
    Species* m_species;

    // breeder of the problem
    Breeder* m_breeder;

    // best individual by far
    Individual m_best_ind;

    // population
    vector<Individual*> m_population;

    // methods to create individuals
    vector<string> m_creator_methods;

    // number of the individuals in the population
    int m_population_size;

    // current generation
    int m_cur_generation;

    // maximum generation
    int m_max_generateions;

    // whether goal is achieved
    bool m_goal_achieved;

    // stat
    clock_t m_start_time;

    // number of generations the best inidividual survived
    int m_generation_survived;

    // random generator
    Random m_random;

    // average fitness in each generation
    //vector<double> m_gen_avg_fitness;

    // best fitness in each generation
    //vector<double> m_gen_best_fitness;
};

} /* namespace ga */
} /* namespace masc */

#endif /* PROBLEM_H_ */
