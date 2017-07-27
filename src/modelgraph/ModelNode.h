//------------------------------------------------------------------------------
//  Copyright 2007-2009 by Jyh-Ming Lien and George Mason University
//  See the file "LICENSE" for more information
//------------------------------------------------------------------------------

#ifndef _MODELGRAPH_NODEL_H_
#define _MODELGRAPH_NODEL_H_

#include <list>
using namespace std;

//////////////////////////////////////////////////////////////////////
// Include ModelGraph Headers.
#include "ModelEdge.h"

class CModelNode  
{
public:
    //////////////////////////////////////////////////////////////////////
    // Constructor/Destructors
    CModelNode(int Key);

    //////////////////////////////////////////////////////////////////////
    // Access Function
    CModelEdge * getEdge(const CModelNode& nb) const
    {
        int Key=nb.m_Key;
        //linear search
        for(list<CModelEdge *>::const_iterator i=m_Edges.begin();i!=m_Edges.end();i++)
            if( (*i)->isEndPt(Key)==true ) return *i;
        return NULL;
    }

    const list<CModelEdge *>& getEdges() const
    {
        return m_Edges;
    }

    void addEdge( CModelEdge * e){ 
        if(e==NULL) return;
        m_Edges.push_back(e);
    }

    int getKey() const { return m_Key; }

    //////////////////////////////////////////////////////////////////////
    // Private Stuff
private:
    
    int m_Key;
    list<CModelEdge*> m_Edges;
};


#endif // _MODELGRAPH_NODEL_H_
