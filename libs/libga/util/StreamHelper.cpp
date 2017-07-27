/*
 * StreamHelper.cpp
 *
 *  Created on: Mar 16, 2015
 *      Author: zhonghua
 */

#include "StreamHelper.h"
#include <cstring>
using namespace std;

namespace masc {
namespace ga {
namespace util {

vector<string> StreamHelper::tokenize(istream& in)
{
    if(in.eof()) return {};

    const int size = 1024;
    char * tmp = new char[size];
    while(!in.eof())
    {
        in.getline(tmp, size);//read lines from a file into strings
        string line = tmp;

        // skip empty line and comments (start with a '#' symbol)
        if(line.length() > 0 && line[0] != '#')
        {
            break;
        }
    }
    auto tok = tokenize(tmp, " =\t[]()<>,");
    delete [] tmp;
    return tok;
}

vector<string> StreamHelper::tokenize(char * line, const char * ignore)
{
    vector<string> tokens;
    char * tok = strtok(line, ignore);
    while (tok != NULL)
    {
        tokens.push_back(tok);
        tok = strtok(NULL, ignore);
    }
    return tokens;
}

} /* namespace until */
} /* namespace ga */
} /* namespace masc */
