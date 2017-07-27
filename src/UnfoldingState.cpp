/*
 * UnfoldingState.cpp
 *
 *  Created on: Nov 16, 2016
 *      Author: zxi
 */

#include "UnfoldingState.h"

namespace masc {
namespace unfolding {

const string UnfoldingState::UNKOWN_PROPERTY = "Unknown";

string UnfoldingState::ToString() const {
  string output = "{";

  bool flag = false;

  for (const auto& kv : this->properties_) {
    if (flag)
      output += ", ";
    output += "\"" + kv.first + "\" : ";
    output += "\"" + kv.second + "\"";
    flag = true;
  }

  output += "}";

  return output;
}

void UnfoldingState::print(ostream& out) const {
  out << this->ToString();
}

} /* namespace unfolding */
} /* namespace masc */

ostream& operator<<(ostream& out,
    const masc::unfolding::UnfoldingState& state) {
  state.print(out);
  return out;
}
