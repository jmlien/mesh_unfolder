/*
 * Species.cpp
 *
 *  Created on: Mar 16, 2015
 *      Author: zhonghua
 */

#include "Species.h"
#include <fstream>

#include "Problem.h"

namespace masc {
namespace ga {

Species::Species(Problem* problem) {
    this->m_genome_size = 1;
    this->m_chunk_size = 1;
    this->m_min_gene = 0.0;
    this->m_max_gene = 1.0;

    this->m_mutation_prob = 0.1;
    this->m_crossover_prob = 0.9;

    this->m_crossover_type = CrossoverType::ONE_POINT;
    this->m_mutation_type = MutationType::RESET;

    this->m_mutation_gaussian_stddev = 1.0;

    this->m_problem = problem;
}

Species::~Species() {
    this->m_problem = nullptr;
}

// setup from tokens
bool Species::setup(const vector<string>& tokens)
{
    for(int i=1;i<tokens.size();++i)
    {
        const string& token = tokens[i];

        if(token == "genome-size")
        {
            this->m_genome_size = stoi(tokens[++i]);
        }
        else if(token == "chunck-size")
        {
            this->m_chunk_size = stoi(tokens[++i]);
        }
        else if(token == "min-gene")
        {
            this->m_min_gene = stof(tokens[++i]);
        }
        else if(token == "max-gene")
        {
            this->m_max_gene = stof(tokens[++i]);
        }
        else if(token == "mutation-prob")
        {
            this->m_mutation_prob = stof(tokens[++i]);
        }
        else if(token == "mutation-type")
        {
            this->m_mutation_type = MutationTypeConverter::parse(tokens[++i]);
        }
        else if(token == "crossover-type")
        {
            this->m_crossover_type = CrossoverTypeConverter::parse(tokens[++i]);
        }
        else if(token == "crossover-prob")
        {
          this->m_crossover_prob = stof(tokens[++i]);
        }
        else if(token == "mutation-stddev")
        {
            this->m_mutation_gaussian_stddev = stof(tokens[++i]);
        }
        else
        {
            cerr<<"- [GA] ! Error in Species::setup line:"<<__LINE__<<". Unknown token = "<<token<<endl;
            return false;
        }
    }

    return true;
}

void Species::print(ostream& out)
{
    char buf[65536];

    sprintf(buf, "species min-gene=%f max-gene=%f genome-size=%d chunck-size=%d mutation-type=%s mutation-prob=%f crossover-type=%s crossover-prob=%f mutation-stddev=%f",
            this->m_min_gene,
            this->m_max_gene,
            this->m_genome_size,
            this->m_chunk_size,
            MutationTypeConverter::toString(this->m_mutation_type).c_str(),
            this->m_mutation_prob,
            CrossoverTypeConverter::toString(this->m_crossover_type).c_str(),
            this->m_crossover_prob,
            this->m_mutation_gaussian_stddev);


    out<<"- [GA] "<<buf<<endl;
}


Individual* Species::createNewIndvidual()
{
    auto ind = new Individual(this->m_genome_size, this);

    return ind;
}

//////////////////////////////////////////
// access
//////////////////////////////////////////

void Species::setMinGene(const float min_gene)
{
    this->m_min_gene = min_gene;
}

const float Species::getMinGene() const
{
    return this->m_min_gene;
}

void Species::setMaxGene(const float max_gene)
{
    this->m_max_gene = max_gene;
}

const float Species::getMaxGene() const
{
    return this->m_max_gene;
}

const int Species::getGenomeSize() const
{
    return this->m_genome_size;
}

void Species::setGenomeSize(int genome_size)
{
    this->m_genome_size = genome_size;
    //cerr<<"!Override! Species::genomeSize = "<<genome_size<<endl;
}

const int Species::getChunkSize() const
{
    return this->m_chunk_size;
}

void Species::setChunkSize(int chunk_size)
{
    this->m_chunk_size = chunk_size;
    //cerr<<"!Override! Species::chunkSize = "<<chunk_size<<endl;
}

void Species::setMutationProb(const float prob)
{
    this->m_mutation_prob = prob;
}

const float Species::getMutationProb() const
{
    return this->m_mutation_prob;
}

void Species::setMutationGaussianStddev(const double stddev)
{
    this->m_mutation_gaussian_stddev = stddev;
}

const double Species::getMutationGaussianStddev() const
{
    return this->m_mutation_gaussian_stddev;
}

void Species::setCrossoverProb(const float prob)
{
    this->m_crossover_prob = prob;
}

const float Species::getCrossoverProb() const
{
    return this->m_crossover_prob;
}

const MutationType Species::getMutationType() const
{
    return this->m_mutation_type;
}

const CrossoverType Species::getCrossoverType() const
{
    return this->m_crossover_type;
}

Problem* Species::getProblem()
{
    return m_problem;
}
} /* namespace ga */
} /* namespace masc */
