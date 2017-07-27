//------------------------------------------------------------------------------
//  Copyright 2007-2009 by Jyh-Ming Lien and George Mason University
//  See the file "LICENSE" for more information
//------------------------------------------------------------------------------

#include "ModelGraph.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////
CModelFacet::CModelFacet()
{
    m_Edge[0]=m_Edge[1]=m_Edge[2]=NULL;
}

CModelGraph::CModelGraph()
{
    this->m_pEdge = NULL;
    this->m_EdgeSize = 0;
    this->m_Tail = NULL;
    this->m_Model = NULL;
}

CModelGraph::~CModelGraph()
{
    //Free edges
    while( m_pEdge!=NULL ){
        CModelEdge * tmp=m_pEdge->getNext();
        delete m_pEdge;
        m_pEdge=tmp;
    }
    //Free Nodes
    m_Nodes.clear();
}

//////////////////////////////////////////////////////////////////////
// init
//////////////////////////////////////////////////////////////////////

//convert from tri to graph
bool CModelGraph::doInit(model* m)
{
    m_Tail=NULL; //tail of m_pEdge

    //create nodes
    m_Nodes.reserve(m->v_size);
    for( uint iP=0;iP<m->v_size;iP++ )
        m_Nodes.push_back(CModelNode(iP));

    //create edge from triangles
    unsigned int eid=0;
    vector<triangle>& tri=m->tris;
    m_Facets.reserve(m->t_size);
    for( uint iT=0; iT<m->t_size; iT++ ){
        CModelFacet F;
        for(int iV=0;iV<3;iV++) F.vid[iV]=tri[iT].v[iV];
        Point3d& v1=m->vertices[F.vid[0]].p;
        Point3d& v2=m->vertices[F.vid[1]].p;
        Point3d& v3=m->vertices[F.vid[2]].p;
        F.n=((v2-v1)%(v3-v1)).normalize();

        for( int iE=0;iE<3;iE++ ){ //for each edge
            int a=tri[iT].v[iE]; int b=tri[iT].v[(iE+1)%3]; //int c=tri[iT][(iE+2)%3];
            CModelEdge * e=m_Nodes[a].getEdge(m_Nodes[b]);
            if( e==NULL ){
                if( (e=new CModelEdge(a,b,eid++))==NULL ) return false;
                m_Nodes[a].addEdge(e); m_Nodes[b].addEdge(e);
                //push new edge to end of list
                if(m_pEdge==NULL) m_pEdge=e;
                else m_Tail->setNext(e);
                m_Tail=e; m_EdgeSize++;
                e->addFacet(iT);
                e->getV()=(m->vertices[b].p-m->vertices[a].p).normalize();
                e->addInNormal((F.n%e->getV()).normalize());
            } //end of if
            else{
                e->addFacet(iT);
                e->addInNormal((F.n%(-e->getV())).normalize());
            }

            F.m_Edge[iE]=e;
        }//end of for(iE)

        m_Facets.push_back(F);

    }//end of for(iT)
    m_Tail->setNext(NULL);
    m_Model=m;

    return true;
}






//////////////////////////////////////////////////////////////////////
// compute the normal of points on 1-d manifold
//////////////////////////////////////////////////////////////////////
//
//void CModelGraph::compute_1_manifold_normal()
//{
//  //for each node i
//  for(vector<CModelNode>::iterator i=m_Nodes.begin();i!=m_Nodes.end();i++){
//      const list<CModelEdge *>& edges=i->getEdges();
//      Vector3d& pn=m_Model->pts[i->getKey()].n;
//      pn.set(0,0,0);
//      for(list<CModelEdge*>::const_iterator j=edges.begin();j!=edges.end();j++){
//          CModelEdge* edge=*j;
//          if(edge->getFacet1()>=0) continue; //not boundary
//          const point& s=m_Model->pts[edge->getStartPt()];
//          const point& e=m_Model->pts[edge->getEndPt()];
//          Vector3d v=e.p-s.p;
//          Vector3d n(-v[1],v[0],0);
//          const Vector3d& fn=m_Facets[edge->getFacet0()].n;
//          if( (n%v)*fn<0 ) n=-n;
//          pn=pn+n.normalize();
//      }
//      pn=pn.normalize();
//  }//end for i
//}
//
////////////////////////////////////////////////////////////////////////
//// split
////////////////////////////////////////////////////////////////////////
//
//
//
////generate points on a 2-manifold
//void CModelGraph::split_1_manifold(float d, vector<point>& ext)
//{
//   CModelEdge* ptr=m_pEdge;
//    int VID=m_Model->pt_size;
//    while(ptr!=NULL){
//        const point& s=m_Model->pts[ptr->getStartPt()];
//        const point& e=m_Model->pts[ptr->getEndPt()];
//  if(ptr->getFacet1()<0)
//          ptr->subdivide(d,s,e,VID,ext);
//        ptr=ptr->getNext();
//    }
//}
//
//
////generate points on a 2-manifold
//void CModelGraph::split_2_manifold
//(float d, vector<point>& ext, int edgescale)
//{
//    CModelEdge* ptr=m_pEdge;
//    int VID=m_Model->pt_size;
//    while(ptr!=NULL){
//        const point& s=m_Model->pts[ptr->getStartPt()];
//        const point& e=m_Model->pts[ptr->getEndPt()];
//      int S=1;
//      if(edgescale!=1){
//          if(ptr->getFacet1()<0){
//              S=edgescale;
//          }
//          else{
//              const Vector3d& n1=m_Facets[ptr->getFacet0()].n;
//              const Vector3d& n2=m_Facets[ptr->getFacet1()].n;
//              if(n1*n2<0.965f)
//                  S=edgescale;
//          }
//      }
//      ptr->subdivide(d/S,s,e,VID,ext);
//        ptr=ptr->getNext();
//    }
//
//    //split all existing E
//  list<CModelFacet> tmp;
//  tmp.insert(tmp.end(),m_Facets.begin(),m_Facets.end());
//    //int index=0;
//    while(!tmp.empty()){
//
//        //index++;
//
//        CModelFacet f=tmp.front();
//        tmp.pop_front();
//
//        bool no_split=false;
//        for(int j=0;j<3;j++){
//            if(f.m_Edge[j]->getSubs().empty())
//                no_split=true;
//        }
//
//        if(!no_split){ //now...create 4 more triangles, 3 more edges
//
//            CModelEdge * e1=f.m_Edge[0]; e1->createKids();
//            CModelEdge * e2=f.m_Edge[1]; e2->createKids();
//            CModelEdge * e3=f.m_Edge[2]; e3->createKids();
//
//            int mid[3]={e1->getMidPt(),e2->getMidPt(),e3->getMidPt()};
//
//          CModelEdge * ne1, * ne2, * ne3;
//          {
//              const point& midpt0=getPoint(mid[0],ext);
//              const point& midpt1=getPoint(mid[1],ext);
//              const point& midpt2=getPoint(mid[2],ext);
//
//              ne1=new CModelEdge(mid[0],mid[2]);
//              ne2=new CModelEdge(mid[1],mid[0]);
//              ne3=new CModelEdge(mid[2],mid[1]);
//
//              vector<point> tmp;
//              ne1->subdivide(d,midpt0,midpt2,VID,tmp);
//              ne2->subdivide(d,midpt1,midpt0,VID,tmp);
//              ne3->subdivide(d,midpt2,midpt1,VID,tmp);
//              ext.insert(ext.end(),tmp.begin(),tmp.end());
//          }
//
//            CModelFacet f1,f2,f3,f4; 
//          {
//              CModelEdge * ke[6]={e1->getKid0(),e1->getKid1(),
//                  e2->getKid0(),e2->getKid1(),
//                  e3->getKid0(),e3->getKid1()};
//
//              if(e1->getStartPt()!=f.vid[0]) swap(ke[0],ke[1]);
//              if(e2->getStartPt()!=f.vid[1]) swap(ke[2],ke[3]);
//              if(e3->getStartPt()!=f.vid[2]) swap(ke[4],ke[5]);
//
//              f1.m_Edge[0]=ke[0]; f1.m_Edge[1]=ne1; f1.m_Edge[2]=ke[5];
//              f1.vid[0]=f.vid[0]; f1.vid[1]=mid[0]; f1.vid[2]=mid[2];
//
//              f2.m_Edge[0]=ke[1]; f2.m_Edge[1]=ke[2]; f2.m_Edge[2]=ne2;
//              f2.vid[0]=mid[0]; f2.vid[1]=f.vid[1]; f2.vid[2]=mid[1];
//
//              f3.m_Edge[0]=ne3; f3.m_Edge[1]=ke[3]; f3.m_Edge[2]=ke[4];
//              f3.vid[0]=mid[2]; f3.vid[1]=mid[1]; f3.vid[2]=f.vid[2];
//
//              f4.m_Edge[0]=ne1; f4.m_Edge[1]=ne2; f4.m_Edge[2]=ne3;
//              f4.vid[0]=mid[0]; f4.vid[1]=mid[1]; f4.vid[2]=mid[2];
//          }
//
//            //--------------------------------------------------------
//            tmp.push_back(f1);
//            tmp.push_back(f2);
//            tmp.push_back(f3);
//            tmp.push_back(f4);
//            m_Tail->setNext(ne1); e1->setNext(ne2);
//            ne2->setNext(ne3);     m_Tail=ne3;
//        }
//        m_Tail->setNext(NULL);
//
//    }//end i
//
//}
//
