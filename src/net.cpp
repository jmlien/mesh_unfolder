#include "net.h"
#include "unfolder.h"

namespace masc {

Net::Net(model * m, const set<uint>& crease_edges, Net * p1, Net * p2)
{
    m_orig=m;
    m_subm=NULL;
    m_parents[0]=p1;
    m_parents[1]=p2;
    m_creases=crease_edges;

    set<uint> faceset;
    for(auto eid:crease_edges)
    {
      edge & e= m->edges[eid];
      faceset.insert(e.fid.begin(), e.fid.end());
    }

    m_faces.insert(m_faces.end(), faceset.begin(), faceset.end());
    build_subm();
    unfold_subm();
}

Net::Net(Net * n)
{
  m_orig=n->m_orig;
  m_subm=n->m_subm;
  m_faces=n->m_faces;
  m_subm_crease_lines=n->m_subm_crease_lines;
  m_creases=n->m_creases;
  m_overlap_pairs=n->m_overlap_pairs;
}

Net::Net(model * m, int fid, Net * p1)
{
    m_orig=m;
    m_subm=NULL;
    m_parents[0]=p1;
    m_parents[1]=NULL;
    m_faces.push_back(fid);
    build_subm();
    unfold_subm();
}

Net::~Net()
{
    //if(m_subm!=NULL) delete m_subm;
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
Net * Net::split(uint eid, uint subm_seed_fid, vector<bool>& face_visited)
{
  set<uint> newcreases;
  list<uint> open; open.push_back(subm_seed_fid);
  face_visited[subm_seed_fid]=true;

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

  if(!newcreases.empty())
    return new Net(m_orig, newcreases, this);
  else //the only face is seed_fid
  {
    return new Net(m_orig, m_subm->tris[subm_seed_fid].source_fid, this);
  }
}

//merge two nets at a given eid
Net * Net::merge(Net * n2, uint eid)
{
  set<uint> merged_creases=this->m_creases;
  merged_creases.insert(n2->m_creases.begin(),n2->m_creases.end());
  merged_creases.insert(eid);
  return new Net(m_orig, merged_creases, this, n2);
}

void Net::build_subm()
{
  m_subm=shared_ptr<model>(m_orig->create_submodel(m_faces));
  assert(m_subm);
}

//unfold m_subm into m_unfolding
void Net::unfold_subm()
{
    //create an unfolder
    Config config;
    config.quite=true;
    config.record_overlap=(m_parents[0]==NULL);
    Unfolder * unfolder=toUnfolder(config, config.record_overlap );

    //remember results
    //this->m_subm_unfolded=unfolder->getUnfoldedMesh();
    const bool * crease=unfolder->getSelectedEdges();
    this->m_subm_crease_lines=vector<bool>(crease,crease+m_subm->e_size);

    if(m_parents[0]==NULL)//no parents, overlaps are created in unfolder
    {
      auto & overlaps = unfolder->getOverlppingFacePairs();
      for(uint a=0;a<m_subm->t_size;a++)
      {
        for(uint b : overlaps[a])
        {
          if(a>=b) continue; //this will be handled by (b,a)
          this->m_overlap_pairs.insert(make_pair(m_subm->tris[a].source_fid,m_subm->tris[b].source_fid));
        }//end for b
      }//end for a
    }
    else //net created by split or merge
    {
       unordered_map<int,uint> fidmap;
       fidmap.reserve(m_subm->t_size);
       for(uint i=0;i<m_subm->t_size;i++)
       {
         fidmap[m_subm->tris[i].source_fid]=i;
       }

       //
       if(m_parents[1]==NULL){ //net created by split
          //
          for(auto& overlap : m_parents[0]->m_overlap_pairs)
          {
            int f1=overlap.first;
            int f2=overlap.second;
            if( fidmap.find(f1)==fidmap.end() ) continue;
            if( fidmap.find(f2)==fidmap.end() ) continue;
            if(unfolder->checkOverlapNew(fidmap[f1],fidmap[f2]))
              this->m_overlap_pairs.insert(overlap);
          }
       }
       else{ //m_parents[0]!=NULL && m_parents[1]!=NULL
          //net created by merging
          this->m_overlap_pairs=m_parents[0]->m_overlap_pairs;
          this->m_overlap_pairs.insert(m_parents[1]->m_overlap_pairs.begin(), m_parents[1]->m_overlap_pairs.end());
          for(auto& f1 : m_parents[0]->m_faces)
          {
            for(auto& f2 : m_parents[1]->m_faces)
            {
              if(unfolder->checkOverlapNew(fidmap[f1],fidmap[f2]))
                this->m_overlap_pairs.insert( (f1<f2)?make_pair(f1,f2):make_pair(f2,f1));
            }//end f2
          }//end f1
       }
    }

    delete unfolder;
}

Unfolder * Net::toUnfolder(const Config & config, bool checkoverlap)
{
  if(this->m_subm==NULL) build_subm();

  //create an unfolder
  Unfolder * unfolder = new Unfolder(this->m_subm.get(),config);

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

  return unfolder;
}

}//end namespace masc
