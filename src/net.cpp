#include "net.h"
#include "unfolder.h"

namespace masc {

Net::Net(model * m, const set<uint>& crease_edges)
{
    m_orig=m;
    m_subm=NULL;
    m_creases=crease_edges;

    set<uint> faceset;
    for(auto eid:crease_edges)
    {
      edge & e= m->edges[eid];
      faceset.insert(e.fid.begin(), e.fid.end());
    }

    m_faces.insert(m_faces.end(), faceset.begin(), faceset.end());
    build_subm();
}

Net::Net(model * m, int fid)
{
    m_orig=m;
    m_subm=NULL;
    m_faces.push_back(fid);
    build_subm();
}

Net::~Net()
{
    if(m_subm!=NULL) delete m_subm;
}

//split the net along the edge eid
pair<Net *, Net*> Net::split(uint eid)
{
  const pair<Net *, Net*> invalid_ans=make_pair(this, (Net*)NULL);
  edge & split_e= this->m_orig->edges[eid];
  if(split_e.fid.size()<2) return invalid_ans;  //nothing to split
  if(split_e.fid.size()>2){
    cerr<<"! Error: Net::split cannot split a non-manifold edge"<<endl;
    return invalid_ans;
  }

  //now we can split faces into two lists
  vector<bool> face_visited(m_subm->t_size,false);
  list<Net*> new_nets;
  for(uint i=0;i<m_subm->t_size;i++)
  {
    if(face_visited[i]==false)
    {
      Net * tmp=split(eid,i,face_visited);
      new_nets.push_back(tmp);
    }
  }
  //well make sure there are only two
  if(new_nets.size()!=2)
  {
    cerr<<"! Error: Net::split: there should only be two new nets but "
        <<new_nets.size()<<" is created"<<endl;
    return invalid_ans;
  }

  return make_pair(new_nets.front(),new_nets.back());
}

//split the net using a seed_fid and cut edge id (eid)
Net * Net::split(uint eid, uint seed_fid, vector<bool>& face_visited)
{
  set<uint> newcreases;
  list<uint> open; open.push_back(seed_fid);
  face_visited[seed_fid]=true;

  while(open.empty()==false)
  {
    uint fid=open.front();
    open.pop_front();
    triangle & tri=m_subm->tris[fid];

    for(uint sub_eid : tri.e) //visit neighbors
    {
      //check if not a crease line....cannot cross
      if( this->m_subm_crease_lines[sub_eid]==false ) continue;
      uint nid=m_subm->edges[sub_eid].otherf(fid);
      if(face_visited[nid]) continue; //visited before
      triangle & nei=m_subm->tris[nid];
      uint orig_eid=m_orig->getEdgeIdByFids(tri.source_fid, nei.source_fid);
      if( orig_eid==eid )
        continue; //this is the edge that we want to cut
      //now we get to the next face
      open.push_back(nid);
      newcreases.insert(orig_eid);
      face_visited[nid]=true;
    }
  }//end while

  cout<<"newcreases size="<<newcreases.size()<<endl;

  if(!newcreases.empty())
    return new Net(m_orig, newcreases);
  else //the only face is seed_fid
    return new Net(m_orig, seed_fid);
}

//merge two nets at a given eid
Net * Net::merge(Net * n2, uint eid)
{
  set<uint> merged_creases=this->m_creases;
  merged_creases.insert(n2->m_creases.begin(),n2->m_creases.end());
  merged_creases.insert(eid);
  return new Net(m_orig, merged_creases);
}

void Net::build_subm()
{
  if(m_subm!=NULL) delete m_subm;

  m_subm=m_orig->create_submodel(m_faces);
  assert(m_subm);
  unfold_subm();
}

//unfold m_subm into m_unfolding
void Net::unfold_subm()
{
    //create an unfolder
    Config config;
    config.quite=true;
    Unfolder * unfolder=toUnfolder(config, true);

    //remember results
    this->m_unfolded=unfolder->getUnfoldedMesh();
    const bool * crease=unfolder->getSelectedEdges();
    this->m_subm_crease_lines=vector<bool>(crease,crease+m_subm->e_size);
    auto & overlaps = unfolder->getOverlppingFacePairs();
    for(uint a=0;a<m_subm->t_size;a++)
    {
      for(uint b : overlaps[a])
      {
        if(a>=b) continue; //this will be handled by (b,a)
        //cout<<"overlaps a="<<a<<" b="<<b<<endl;
        this->m_overlap_pairs.insert(make_pair(m_subm->tris[a].source_fid,m_subm->tris[b].source_fid));
      }//end for b
    }//end for a

    delete unfolder;
}

Unfolder * Net::toUnfolder(const Config & config, bool checkoverlap)
{
  //create an unfolder
  Unfolder * unfolder = new Unfolder(this->m_subm,config);

  //create weights so all weights are 1 initiall and 0 for crease_style
  vector<float> weights(m_subm->e_size, 1);
  //for (uint eid : this->m_creases) weights[eid] = 0.0;
  for(uint i=0;i<m_subm->e_size;i++)
  {
    edge & e=m_subm->edges[i];
    if(e.fid.size()!=2) continue;
    triangle& t1=m_subm->tris[e.fid.front()];
    triangle& t2=m_subm->tris[e.fid.back()];
    uint eid=m_orig->getEdgeIdByFids(t1.source_fid,t2.source_fid);
    if(this->m_creases.find(eid)!=this->m_creases.end())
      weights[i] = 0.0;
  }

  //unfold using weights
  unfolder->buildFromWeights(weights, checkoverlap);

  //unfold using weights
  //unfolder.buildFromWeights(weights, true);
  //unfolder.rebuildModel();
  //unfolder.unfoldTo(0.0);
  //unfolder.dumpObj("debug.obj");
  //unfolder.dumpSVG("debug.svg", ExportSVGType::BASIC);

  return unfolder;
}

}//end namespace masc
