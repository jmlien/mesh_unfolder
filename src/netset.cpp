#include "netset.h"
#include <cfloat>

//
//a set of nets spanning the given mesh
//
namespace masc {

NetSet::NetSet(Net * n)
{
  this->m_nets.push_back(n);
  this->m_orig=n->getOriginalModel();
  const auto& net_faces=n->getFaces();
  assert(this->m_orig->t_size==net_faces.size());
  m_fid2net=vector<Net *>(net_faces.size(),n);
}

//copy constructor
NetSet::NetSet(NetSet * ns)
{
  m_orig=ns->m_orig;
  m_nets=ns->m_nets;
  //m_created_nets=ns->m_created_nets;
  m_fid2net=ns->m_fid2net;
  //unordered_map<Net*,Net*> nmap;
  // for(Net * n : ns->m_nets)
  // {
  //   Net * newn=new Net(n);
  //   m_nets.push_back(newn);
  //   m_created_nets.push_back(newn);
  //   nmap[n]=newn;
  // }
  //
  // m_fid2net.resize(ns->m_fid2net.size());
  // for(Net * n : ns->m_fid2net) m_fid2net.push_back(nmap[n]);
}

NetSet::~NetSet()
{
  for(auto net: m_created_nets) delete net;
  m_nets.clear();
  m_created_nets.clear();
}

bool NetSet::split(uint eid) //split a net at a given edge, this edge must be a crease edge
{
  //make sure that eid is a crease edge
  edge & e=m_orig->edges[eid];
  if(e.type=='b') return false;
  assert(e.fid.size()>1);
  uint fid1=e.fid.front();
  uint fid2=e.fid.back();

  Net * net1=m_fid2net[fid1];
  Net * net2=m_fid2net[fid2];

  if(net1!=net2) return false; //cut line
  if( net1->is_crease(eid)==false ) return false; //not a crease

  auto new_nets=net1->split(eid);
  if(new_nets.second==NULL) return false; //sometime wrong happended

  //remove net1 & net2 and add net3
  replace(new_nets.first);
  replace(new_nets.second);
  m_nets.remove(net1);
  m_nets.push_back(new_nets.first);
  m_nets.push_back(new_nets.second);

  m_created_nets.push_back(new_nets.first); //for garbage collection
  m_created_nets.push_back(new_nets.second);//for garbage collection

  return true;
}

bool NetSet::merge(uint eid) //merge two nets at a given edge, this edge must be a cut edge
{
  //make sure that eid is a cut edge between two nets
  edge & e=m_orig->edges[eid];
  if(e.type=='b') return false;
  assert(e.fid.size()>1);
  uint fid1=e.fid.front();
  uint fid2=e.fid.back();
  Net * net1=m_fid2net[fid1];
  Net * net2=m_fid2net[fid2];
  if(net1==net2){
    cerr<<"! Error: NetSet::merge: Trying merge the same net ("<<*net1<<") at edge "<<eid<<endl;
    return false; //same net...
  }

  Net * net3=m_fid2net[fid1]->merge(m_fid2net[fid2], eid);

  //remove net1 & net2 and add net3
  replace(net3);
  m_nets.remove(net1);
  m_nets.remove(net2);
  m_nets.push_back(net3);
  m_created_nets.push_back(net3); //for garbage collection

  return true;
}

void NetSet::sharedCutEdges(Net * n1, Net * n2, set<uint>& shared)
{
  //go through edges of n1
  for(uint f1 : n1->getFaces())
  {
    const triangle & tri = m_orig->tris[f1];
    for(int i=0;i<3;i++)
    {
      uint eid=tri.e[i];
      if(m_orig->edges[eid].type=='b') continue; //border edge
      uint f2 = m_orig->edges[eid].otherf(f1);
      if(f2==((uint)-1)) continue; //invalid f2....
      if( m_fid2net[f2]!=n2 ) continue; //not shared between n1 and n2
      shared.insert(eid);
    }//end for i
  }//end for f1

  //go through edges of n2
  for(uint f1 : n2->getFaces())
  {
    const triangle & tri = m_orig->tris[f1];
    for(int i=0;i<3;i++)
    {
      uint eid=tri.e[i];
      if(m_orig->edges[eid].type=='b') continue; //border edge
      uint f2 = m_orig->edges[eid].otherf(f1);
      if(f2==((uint)-1)) continue; //invalid f2....
      if( m_fid2net[f2]!=n1 ) continue; //not shared between n1 and n2
      shared.insert(eid);
    }//end for i
  }//end for f1
}

// evaluate this netset
// return the evaluation. The higher value is worse. 0 is the best.
double NetSet::evaluate()
{
  double sum=0;
  for(Net * net: m_nets) sum+=net->getOverlaps().size();
  return sum;
}

// evaluate merging of two nets at a given edge, this edge must be a cut edge
// return the evaluation. The higher value is worse. 0 is the best.
double NetSet::evaluate_merge(uint eid)
{
  //make sure that eid is a cut edge between two nets
  edge & e=m_orig->edges[eid];
  if(e.type=='b') return DBL_MAX; //can't be merged

  assert(e.fid.size()>1);
  uint fid1=e.fid.front();
  uint fid2=e.fid.back();
  Net * net1=m_fid2net[fid1];
  Net * net2=m_fid2net[fid2];
  if(net1==net2) return false; //same net...

  Net * net3=m_fid2net[fid1]->merge(m_fid2net[fid2], eid);

  //the score can be improved to consider penetration depth
  double score = net3->getOverlaps().size()*1.0;

  delete net3;

  return score;
}

//go through faces fid in net and replace m_fid2net[fid] with net
void NetSet::replace(Net * net)
{
  for(uint fid: net->getFaces())
  {
    m_fid2net[fid]=net;
  }
}

}//namespace masc
