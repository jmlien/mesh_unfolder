/*
 * CD.cpp
 *
 *  Created on: Apr 26, 2015
 *      Author: zhonghua
 */

#include "CD.h"
#include "unfolder.h"
#include "RAPID/RAPID.h"

namespace masc {

CD::CD(Unfolder* unfolder) {
    this->m_unfodler = unfolder;
}

CD::~CD() {
    // do not delete the ptr...
    this->m_unfodler = nullptr;
}

//////////////////////////////////////////////

RAPID_CD::RAPID_CD(Unfolder* unfolder) : CD(unfolder) {

}

// do nothing
RAPID_CD::~RAPID_CD() {}


int RAPID_CD::hasCollision(bool first_contact) {

    const auto& mesh = this->m_unfodler->getUnfolded();
    auto model = this->m_unfodler->getModel();

    RAPIDres res;

    int count = 0;

    // rotation
    double R[3][3] = {
            {1,0,0},
            {0,1,0},
            {0,0,1}
    };

    // translation
    double T[3] = {0,0,0};

    const auto FLAG = first_contact ? RAPID_FIRST_CONTACT : RAPID_ALL_CONTACTS;

    // build the body model
    auto body = this->buildBody(-1);

    for(int i=0;i<mesh.size();++i)
    {
        // build a triangle model
        auto tri = this->buildTriangle(i);

        RAPID_Collide(res,
                R,T,body,
                R,T,tri,
                FLAG);

        bool self_intersected = false;

        if(res.RAPID_num_contacts>0)
        {
            for(auto c=0;c<res.RAPID_num_contacts;++c)
            {
                if (res.RAPID_contact[c].id1 == i)
                    self_intersected = true;
                //else
                //    cout<<" - "<<i<<" <-> "<<res.RAPID_contact[c].id1<<endl;
            }
        }

        if(self_intersected) {
            --res.RAPID_num_contacts;
        };

        count += res.RAPID_num_contacts;

        if(res.RAPID_num_contacts > 0) {
            //cout<<"face "<<i<<" intersected with "<<res.RAPID_num_contacts<<" faces"<<endl;
            model->tris[i].overlapped = true;
        }

        delete tri;

        if(count > 0 && first_contact)
            break;
    }

    delete body;

    return count;
}

// build rapid model for entire unfolding...
RAPID_model* RAPID_CD::buildBody(const int fid)
{
    const auto& mesh = this->m_unfodler->getUnfolded();

    auto model = new RAPID_model();

    model->BeginModel();

    int id = 0;
    for(const auto& f : mesh)
    {
        if(id != fid)
            model->AddTri(f[0].second.get(), f[1].second.get(), f[2].second.get(), id);

        ++id;
    }

    model->EndModel();

    return model;
}

// build rapid model for a triangle
RAPID_model* RAPID_CD::buildTriangle(const int fid)
{
    auto model = new RAPID_model();

    const auto& mesh = this->m_unfodler->getUnfolded();
    const auto& f = mesh[fid];

    model->BeginModel();

    model->AddTri(f[0].second.get(), f[1].second.get(), f[2].second.get(), 0);

    model->EndModel();

    return model;
}

} /* namespace masc */
