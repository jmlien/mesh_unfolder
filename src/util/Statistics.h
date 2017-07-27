/*
 * Statistics.h
 *
 *  Created on: Oct 25, 2016
 *      Author: zxi
 */

#ifndef SRC_UTIL_STATISTICS_H_
#define SRC_UTIL_STATISTICS_H_

#include <vector>
#include <limits>
using namespace std;

namespace masc {
namespace util {

template<typename T, typename = typename std::enable_if<
    std::is_arithmetic<T>::value, T>::type>
class Statistics {
public:
  Statistics() {
    min_ = std::numeric_limits<T>::max();
    max_ = std::numeric_limits<T>::min();

    median_ = T(0);
    mean_ = 0.0;
    stddev_ = 0.0;
    dirty_ = false;

    sum_ = 0;
  }

  Statistics(const vector<T>& data) :
      Statistics() {
    for (const auto e : data) {
      this->Add(e);
    }
  }

  ~Statistics() {
  }

  void Reset() {
    data_.clear();
    min_ = std::numeric_limits<T>::max();
    max_ = std::numeric_limits<T>::min();

    median_ = T(0);
    mean_ = 0.0;
    stddev_ = 0.0;
    sum_ = 0.0;
    dirty_ = false;
  }

  const unsigned int Size() {
    return data_.size();
  }

  // Add a new value.
  // Min, Max, mean will be updated.
  void Add(const T& val) {
    data_.push_back(val);
    min_ = std::min(min_, val);
    max_ = std::max(max_, val);
    sum_ += val;
    mean_ = sum_ / data_.size();
    dirty_ = true;
  }

  // Compute the median
  void Stat() {

    if (data_.size() == 0) {
      return;
    }

    std::sort(data_.begin(), data_.end());
    median_ = data_[data_.size() / 2];
    dirty_ = false;

    double varience = 0.0;
    for (const T& v : data_) {
      varience += (v - mean_) * (v - mean_);
    }

    if (data_.size() > 0) {
      varience /= data_.size();
      stddev_ = std::sqrt(varience);
    }
  }

  const T& Min() const {
    return min_;
  }

  const T& Max() const {
    return max_;
  }

  const double Mean() const {
    return mean_;
  }

  const T& Sum() const {
    return sum_;
  }

  const T& Median() const {
    return median_;
  }

  const double StandardDeviation() const {
    return stddev_;
  }

private:
  vector<T> data_;
  T min_;
  T max_;
  T median_;
  double sum_;
  double mean_;
  double stddev_;
  bool dirty_;
};

template<typename T, typename = typename std::enable_if<
    std::is_arithmetic<T>::value, T>::type>
ostream& operator<<(ostream& out, Statistics<T>& s) {
  out << s.Size() << " " << s.Min() << " " << s.Mean() << " " << s.Max() << " "
      << s.Sum() << " " << s.Median() << " " << s.StandardDeviation() << endl;
  return out;
}

} /* namespace util */
} /* namespace masc */

#endif /* SRC_UTIL_STATISTICS_H_ */
