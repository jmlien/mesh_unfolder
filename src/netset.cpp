#include "netset.h"

//
//a set of nets spanning the given mesh
//
namespace masc {

NetSet::NetSet(Net * n)
{
  this->m_nets.push_back(n);
  this->m_orig=n->getOriginalModel();

  const auto& net_faces=n->getFaces();
  m_fid2net=vector<Net *>(net_faces.size(),n);
}

NetSet::~NetSet()
{
  for(auto net: m_nets) delete net;
  m_nets.clear();
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
  if(net1==net2) return false; //same net...

  Net * net3=m_fid2net[fid1]->merge(m_fid2net[fid2], eid);

  //remove net1 & net2 and add net3
  replace(net3);
  m_nets.remove(net1);
  m_nets.remove(net2);
  m_nets.push_back(net3);

  return true;
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
