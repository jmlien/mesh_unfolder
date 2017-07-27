/*
 * UnfoldingState.h
 *
 *  Created on: Nov 16, 2016
 *      Author: zxi
 */

#ifndef SRC_UNFOLDINGSTATE_H_
#define SRC_UNFOLDINGSTATE_H_

#include <vector>
#include <unordered_map>
#include <string>
#include <iostream>
using namespace std;

namespace masc {
namespace unfolding {

// Represent an unfolding state (could have overlapping)
class UnfoldingState {
public:
  UnfoldingState() {
  }
  UnfoldingState(const vector<float>& weights) :
      weights_(weights) {
  }

  void AddProperty(const string& key, const string& value) {
    properties_[key] = value;
  }

  const vector<float>& GetWeights() const {
    return weights_;
  }

  const string& GetProperties(const string& key) const {
    if (properties_.count(key)) {
      return properties_.at(key);
    }
    return UNKOWN_PROPERTY;
  }

  // Get a string representation of the properties;
  string ToString() const;

  virtual void print (ostream& out) const;

  virtual ~UnfoldingState() {
  }

private:
  vector<float> weights_;
  unordered_map<string, string> properties_;
  static const string UNKOWN_PROPERTY;

};

} /* namespace unfolding */
} /* namespace masc */

ostream& operator<<(ostream& out, const masc::unfolding::UnfoldingState& state);

#endif /* SRC_UNFOLDINGSTATE_H_ */
