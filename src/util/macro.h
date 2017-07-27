/*
 * macro.h
 *
 *  Created on: Dec 28, 2016
 *      Author: zxi
 */

#ifndef SRC_UTIL_MACRO_H_
#define SRC_UTIL_MACRO_H_

#define RETURN_IF_FALSE(expr) bool r=(expr); if (!(r)) return (r);
#define CHECK(r) if(!(r)) { cerr << "Check failed!" << endl; exit(-1); }
#define CHECK_EQ(e, a) \
    auto te=(e); \
    auto ta=(a); \
    if(te != ta) { \
      cerr << "Check failed! " << te << "!=" << ta << " in " << __func__ <<  " at "  << __FILE__  << ":" << __LINE__ << endl; \
      exit(-1); \
    } \

#endif /* SRC_UTIL_MACRO_H_ */
