/*
 * Random.h
 *
 *  Created on: Mar 19, 2015
 *      Author: zhonghua
 */

#ifndef RANDOM_H_
#define RANDOM_H_

#include <random>
using namespace std;

namespace masc {
namespace ga {
namespace util {

class Random {
public:
    Random();
    virtual ~Random();
    void setSeed(int seed);
    double nextDoubleUniform();
    double nextDoubleGuasian(double miu=0.0f, double std_dev = 1.0f);
private:
    mt19937 m_uniform_generator;
    mt19937 m_gaussian_generator;

    uniform_real_distribution<double> m_uniform_distribution;
    normal_distribution<double> m_gaussian_distribution;
};

} /* namespace util */
} /* namespace ga */
} /* namespace masc */

#endif /* RANDOM_H_ */
