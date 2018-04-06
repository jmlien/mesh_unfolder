/*
 * Random.cpp
 *
 *  Created on: Mar 19, 2015
 *      Author: zhonghua
 */

#include "Random.h"

#include <iostream>
#include <ctime>

namespace masc {
namespace ga {
namespace util {

Random::Random() {
    // uniform distribution [0,1)
    this->m_uniform_distribution = std::uniform_real_distribution<double>(0, 1);

    // gaussian distribution mean = 0, standard deviation = 1.0
    this->m_gaussian_distribution = std::normal_distribution<double>(0, 1);
}

Random::~Random() {

}


void Random::setSeed(int seed)
{
    cout<<"- [GA] Random::setSeed: seed set to "<<seed<<endl;
    this->m_uniform_generator.seed(seed);
    this->m_gaussian_generator.seed(seed);
    this->m_uniform_distribution.reset();
    this->m_gaussian_distribution.reset();
}

double Random::nextDoubleUniform()
{
    return this->m_uniform_distribution(this->m_uniform_generator);
}

double Random::nextDoubleGuasian(double miu, double std_dev)
{
    return this->m_gaussian_distribution(this->m_gaussian_generator) * std_dev + miu;
}

} /* namespace util */
} /* namespace ga */
} /* namespace masc */
