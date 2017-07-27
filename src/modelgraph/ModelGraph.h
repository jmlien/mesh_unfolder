//------------------------------------------------------------------------------
//  Copyright 2007-2009 by Jyh-Ming Lien and George Mason University
//  See the file "LICENSE" for more information
//------------------------------------------------------------------------------

#ifndef _MODEL_GRAPH_H_
#define _MODEL_GRAPH_H_

#include "ModelNode.h"
#include "ModelEdge.h"
#include "model.h"


struct CModelFacet
{
    CModelFacet();
    int vid[3];
    Vector3d n;
    CModelEdge * m_Edge[3];
};

class CModelGraph  
{
public:

    //////////////////////////////////////////////////////////////////////
    // Constructor/Destructors
    CModelGraph();
    ~CModelGraph();

    //////////////////////////////////////////////////////////////////////
    // Core Function
    //convert from tri to graph
    bool doInit(model* m);

    //////////////////////////////////////////////////////////////////////
    // Access Function

    //edges
    int getEdgeSize() const { return m_EdgeSize; }
    CModelEdge * getEdges(){ return m_pEdge; }

    //facets
    vector<CModelFacet>& getFacets(){ return m_Facets; }

    //////////////////////////////////////////////////////////////////////
    // Private stuff
private:
    int m_EdgeSize;
    CModelEdge * m_pEdge;           //a list of edges
    CModelEdge * m_Tail;            //point to the end
    vector<CModelNode>  m_Nodes;    //an array of nodes
    vector<CModelFacet> m_Facets;   //an array of facets
    model* m_Model;
};


#endif
