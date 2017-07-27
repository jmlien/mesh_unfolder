/*
 * StreamHelper.h
 *
 *  Created on: Mar 16, 2015
 *      Author: zhonghua
 */

#ifndef STREAMHELPER_H_
#define STREAMHELPER_H_

#include <iostream>
#include <string>
#include <vector>
using namespace std;

namespace masc {
namespace ga   {
namespace util {

class StreamHelper {
public:

    // tokenize a line from the given input stream.
    // note: empty line and comments will be skipped.
    static vector<string> tokenize(istream& in);

private:

    // internal implementation of tokenization
    static vector<string> tokenize(char * line, const char * ignore);
};

} /* namespace util */
} /* namespace ga */
} /* namespace masc */

#endif /* STREAMHELPER_H_ */
