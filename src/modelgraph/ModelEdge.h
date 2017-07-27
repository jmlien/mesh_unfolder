//------------------------------------------------------------------------------
//  Copyright 2007-2009 by Jyh-Ming Lien and George Mason University
//  See the file "LICENSE" for more information
//------------------------------------------------------------------------------

#ifndef _MODEL_GRAPH_EDGE_H_
#define _MODEL_GRAPH_EDGE_H_

#include "mathtool/Vector.h"
using namespace mathtool;

#include <vector>
using namespace std;

class CModelEdge  
{
public:
    //////////////////////////////////////////////////////////////////////
    // Constructor/Destructors
    CModelEdge(int start, int end, int eid);

    //////////////////////////////////////////////////////////////////////
    // Access Function
    bool isEndPt( int key ) const{
        for( int i=0;i<2;i++)
            if( key==m_Key[i] ) return true;
        return false;
    }
    int getStartPt()  const { return m_Key[0]; }
    int getEndPt()    const { return m_Key[1]; }

    //List Access
    CModelEdge * getNext() const { return m_Next; }
    void setNext(CModelEdge * e) { m_Next=e; }

    void addFacet(int id){ m_Fid.push_back(id); }

    const vector<int>& getFacets() const {return m_Fid;}

    void addInNormal(const Vector3d& n){ m_InNormal.push_back(n); }
    Vector3d& getInNormal(int id){ return m_InNormal[id]; }
    Vector3d& getV(){ return m_V; }

    int getID() const { return m_id; }
    //////////////////////////////////////////////////////////////////////
    // Private Stuff
private:

    int m_Key[2]; ///< 0->start, 1->end
    vector<int> m_Fid; //facet ids

    Vector3d m_V;  //points from start to end
    vector<Vector3d> m_InNormal; //inface normal

    //list link
    CModelEdge * m_Next;

    int m_id;
};


#endif // _MODEL_GRAPH_EDGE_H_
