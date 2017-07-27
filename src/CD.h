/*
 * CD.h
 *
 *  Created on: Apr 26, 2015
 *      Author: zhonghua
 */

#ifndef CD_H_
#define CD_H_

#include <vector>
using namespace std;

// forward declaration
class RAPID_model;
class Unfolder;


namespace masc {

class CD {
public:
    CD(Unfolder* unfolder);
    virtual ~CD();
    // return number of collisions
    virtual int hasCollision(bool first_contact=false)=0;
protected:
    Unfolder* m_unfodler;
};

class RAPID_CD : public CD {
public:
    RAPID_CD(Unfolder* unfolder);
    virtual ~RAPID_CD();

    // check collision, return number of collisions
    // first_contact : stop at first contact
    virtual int hasCollision(bool first_contact=false) override;
protected:
    // build rapid model for entire unfolding...
    // without fid...
    virtual RAPID_model* buildBody(const int fid);
    // build rapid model for a triangle
    virtual RAPID_model* buildTriangle(const int fid);
};

} /* namespace masc */

#endif /* CD_H_ */
