#include "netsurgent.h"
#include "unfolder.h"
#include <unordered_map>
#include <deque>
#include <glpk.h> //linear programming solver
#include "util/DisjointSets.h"
#include "polygon/bbox2d.h"

namespace masc {

//
//
// NetSurgent
//
//

bool NetSurgent::operate(Unfolder* unfolder, vector<Unfolder*> & operated)
{
  //get config
  const Config & config = unfolder->getConfig();
  //get model and creases to create a net
  //cout<<"create net"<<endl;
  model * m=unfolder->getModel();
  const set<uint>& creases=unfolder->getFoldEdges();
  Net * net=new Net(m,creases);

  //create a net surgent and apply
  //cout<<"create surgent"<<endl;
  NetSurgent * surgent=NULL;
  switch(config.surgery_method)
  {
    case NetSurgery::SET_COVER_SURGERY: surgent=new SetCoverNetSurgent(); break;
    case NetSurgery::TOPOLOGICAL_SURGERY: surgent=new TopologicalNetSurgent(10); break;
    case NetSurgery::BOXING_SURGERY: surgent=new BoxingNetSurgent(config.svg_boxing_width,config.svg_boxing_height); break;
    case NetSurgery::SUBDIVID_SURGERY:
    case NetSurgery::CAGING_SURGERY:
    default:
      cerr<<"! Error: Net Surgery type unsupported"<<endl;
  }

  if(surgent==NULL) return false;

  //cout<<"surgent operates the net"<<endl;
  NetSet * netset = surgent->apply(net);
  if(netset==NULL) return false;

  list<Net*> & nets=netset->getNets();

  //convert nets to a list of new unfolders
  //cout<<"convert nets to unfolders"<<endl;
  int cid=1;
  Config _config=config;
   _config.quite=true;
   _config.baseface=-1;
   _config.random_baseface = false;
  // _config.output_file=config.output_file;
  // _config.filename=config.filename;
  // _config.texture_path=config.texture_path;

  for(Net* net: nets)
  {
    Unfolder * tmp=net->toUnfolder(_config, true);
    tmp->setClusterId(cid++);
    operated.push_back(tmp);
  }

//cout<<"net surgery done"<<endl;

  //done
  return true;
}

//
//
// SetCoverNetSurgent
//
//

SetCoverNetSurgent::SetCoverNetSurgent()
{

}
//
// //convert the net into a net set
// NetSet * SetCoverNetSurgent::apply2(Net * net)
// {
//     //find all all_overlaps
//     set< pair<uint,uint> >  overlaps = net->getOverlaps();
//
//     cout<<"there are "<<overlaps.size()<<" overlaps"<<endl;
//
//     //find shortest paths (sequence of edges) connecting pairs of overlapping Triangles
//     unordered_map<uint, list<pair<uint,uint> > > e2tmap; //map< eid, pair<uint, uint>=<tid1, tid2> >
//     for(const pair<uint,uint>& overlap : overlaps)
//     {
//       list<uint> epath;
//       //cout<<"find path for "<<overlap.first<<","<<overlap.second<<": ";
//       //find path between overlap.first and overlap.second
//       this->find_path(net, overlap.first, overlap.second, epath);
//
//       for(uint e : epath)
//       {
//         e2tmap[e].push_back(overlap);
//       //  cout<<e<<", ";
//       }
//       // cout<<endl;
//     }
//
//     cout<<"finished finding all paths, e2tmap has size="<<e2tmap.size()<<endl;
//
//     //this might be slow
//     //for each edge, we maintain a list of overlapping pair
//     //greedily, we start from an edge with the largest number
//     //of un-resolved pairs
//     cout<<"solve set-cover problem"<<endl;
//     list<uint> new_cuts;
//     while(overlaps.empty()==false)
//     {
//       //find an edge with the largest cover
//       uint best_e;
//       uint best_size=0;
//       for(auto& tmp : e2tmap) //loop through each edge
//       {
//           uint eid = tmp.first;
//           const list<pair<uint,uint> >& conflicts = tmp.second;
//           //count the number of resolved pairs
//           int count=0;
//           for(const pair<uint,uint> & overlap : conflicts)
//           {
//             if(overlaps.find(overlap)!=overlaps.end()) count++;
//           }
//           //
//           if(count>best_size)
//           {
//             best_size=count;
//             best_e=eid;
//           }
//           else if(count==best_size) //resolve the same number of pairs
//           {
//             //pick one with shorter length
//             model * m = net->getOriginalModel();
//             edge & e1 = m->edges[best_e];
//             edge & e2 = m->edges[eid];
//             if(e2.length<e1.length)
//             {
//               best_e=eid;
//             }
//           }
//       }//end for tmp
//
//       //cout<<"best_e="<<best_e<<" best_size="<<best_size<<endl;
//
//       //found the best edge, remember it and remove all resolved pairs from overlaps
//       //cout<<"remove resolved overlaps"<<endl;
//       new_cuts.push_back(best_e);
//       for(const pair<uint,uint> & overlap : e2tmap[best_e])
//       {
//         overlaps.erase(overlap); //these overlaps have been resolved
//       }//end for overlap
//
//
//       //cout<<"remove best_e"<<endl;
//       e2tmap.erase(best_e); //remove best_e from the map
//
//       //cout<<"overlaps size="<<overlaps.size()<<endl;
//     }//end while
//
//     //cut the edges
//     //cout<<"cut nets"<<endl;
//     NetSet * netset = new NetSet(net);
//     for(uint eid : new_cuts)
//     {
//       //cout<<"cut edge "<<eid<<endl;
//       bool r = netset->split(eid);
//       if(!r) cerr<<"! Error: failed to cut edge "<<eid<<endl;
//       //break;
//     }
//
//     //cout<<"surgery done"<<endl;
//     return netset;
// }


//convert the net into a net set
NetSet * SetCoverNetSurgent::apply(Net * net)
{
    //find all all_overlaps
    set< pair<uint,uint> >  overlaps = net->getOverlaps();

    if(overlaps.size()==0){ //no overlaps, Done
      cerr<<"! Warning: SetCoverNetSurgent::apply: net is already flattened."<<endl;
      return NULL;
    }

    //cout<<"there are "<<overlaps.size()<<" overlaps"<<endl;

    //find shortest paths (sequence of edges) connecting pairs of overlapping Triangles
    list< list<uint> > epaths;
    for(const pair<uint,uint>& overlap : overlaps)
    {
      epaths.push_back(list<uint>());
      list<uint> & epath=epaths.back();
      this->find_path(net, overlap.first, overlap.second, epath);
    }

    list<uint> new_cuts;
    SolveSetCover(epaths, new_cuts);

    //cut the edges
    //cout<<"cut nets"<<endl;
    NetSet * netset = new NetSet(net);
    for(uint eid : new_cuts)
    {
      //cout<<"cut edge "<<eid<<endl;
      bool r = netset->split(eid);
      if(!r) cerr<<"! Error: failed to cut edge "<<eid<<endl;
      //break;
    }

    //cout<<"surgery done"<<endl;
    return netset;
}

//find a shortest path connecting f1 and f2
//list<uint>& path contains a list of edges that the path passes through
void SetCoverNetSurgent::find_path(Net * net, uint f1, uint f2, list<uint>& path)
{

  model * m = net->getOriginalModel();
  std::deque<uint> open;

  bool * visited = new bool[m->t_size];
  memset(visited, 0, sizeof(bool)* m->t_size);

  //add root (f1)
  visited[f1] = true;
  open.push_back(f1);
  m->tris[f1].parent_id = f1;

  //loop until all faces in open are handled
  while (open.empty() == false)
  {
    uint fid = open.front();
    open.pop_front();

    if (fid == f2) //path found!
    {
      break; //done
    }

    for (int i=0;i<3;i++)
    {
      uint eid=m->tris[fid].e[i];
      //make sure eid is a crease edge
      if( net->is_crease(eid)==false ) continue;
      uint kid=m->edges[eid].otherf(fid);
      if (visited[kid]) continue;
      open.push_back(kid);
      visited[kid] = true;
      m->tris[kid].parent_id = fid;
    } //end for kid

  } //end while

  int now = f2;
  while (now != f1)
  {
    int next=m->tris[now].parent_id;
    path.push_front( m->getEdgeIdByFids(now,next) );
    now = next;
  }

  delete [] visited;

}


bool SetCoverNetSurgent::SolveSetCover
(const list< list<uint> > & epaths, list<uint>& cuts)
{
    vector<LP_constraints> constaints;
    unordered_map<uint,uint> eid2index;
    vector<int> index2eid;

    glp_term_out(GLP_OFF); //make glpk silient

    //build mapping between edge ids and index (order of edges in LP)
    uint index=0;
    for(const list<uint> & epath : epaths)
    {
      for(uint eid : epath)
      {
        if(eid2index.find(eid)==eid2index.end()){
          eid2index[eid]=index++;
          index2eid.push_back(eid);
        }
      }//end id
    }//end epath

    //build constraints
    for(const list<uint> & epath : epaths)
    {
      LP_constraints c;
      //cout<<"constaint : ";
      for(uint eid : epath){
        //cout<<eid2index[eid]<<" ";
        c.eids.push_back(eid2index[eid]);
      }
      //cout<<endl;
      c.type = GLP_LO;
      c.lower_bound = 1;
      constaints.push_back(c);
    }

    int constraint_size = constaints.size();
    int variable_size = index2eid.size();

    glp_prob * lp = glp_create_prob();
    assert(lp);
    glp_set_prob_name(lp, "lp");
    glp_set_obj_dir(lp, GLP_MIN);
    glp_add_rows(lp, constraint_size);
    glp_add_cols(lp, variable_size);

    //init rows
    int iaja_size = 0;
    int row_id = 1;

    char tmp[64];
    for (auto & c : constaints)
    {
      sprintf(tmp, "r%08d", row_id);
      glp_set_row_name(lp, row_id, tmp);
      glp_set_row_bnds(lp, row_id, c.type, c.lower_bound, c.upper_bound);
      iaja_size += c.eids.size();
      row_id++;
    }

    //init cols
    for (int i = 1; i <= variable_size; i++)
    {
      char tmp[64];
      sprintf(tmp, "s%08d", i);
      glp_set_col_name(lp, i, tmp);
      //MIP
      glp_set_col_kind(lp, i, GLP_BV);
      //MIP
      {
        //use edge length, minimize cut length
        // auto & e = m->edges[i - 1];
        // const auto  & vec = m->vertices[e.vid[0]].p - m->vertices[e.vid[1]].p;
        // glp_set_obj_coef(lp, i, vec.norm());

        //use edge count, minimize the # of components
        glp_set_obj_coef(lp, i, 1);
      }

    }//end i

    //init ia, ja, and ar
    int * ia = new int[1 + iaja_size];
    int * ja = new int[1 + iaja_size];
    double * ar = new double[1 + iaja_size];
    assert(ia && ja && ar);

    int ia_id = 1;
    row_id = 1;
    for (auto & c : constaints)
    {
      for (auto eid : c.eids)
      {
        ia[ia_id] = row_id;
        ja[ia_id] = eid + 1;
        ar[ia_id] = 1;
        ia_id++;
      }//end for j

      row_id++;
    }//end for i


    glp_load_matrix(lp, iaja_size, ia, ja, ar);

    //assert(glp_simplex(lp, NULL) == 0);
    //assert(glp_get_status(lp) == GLP_OPT);

    glp_iocp parm;
    glp_init_iocp(&parm);
    parm.pp_tech = GLP_PP_ALL;
    parm.presolve = GLP_ON;
    parm.clq_cuts = GLP_ON;
    parm.binarize = GLP_ON;
    //parm.cb_func = callback;
    parm.cb_info = this;

    parm.tm_lim = 600000; //600 sec
    int err = glp_intopt(lp, &parm);
    //cout << "err=" << err << endl;

    double z = glp_mip_obj_val(lp);
    //cout << "objective value=" << z << endl;

    //get mip status
    int glp_prim_stat = glp_mip_status(lp);

    switch (glp_prim_stat)
    {
      case GLP_OPT: cout << "solution is optimal;" << endl; break;
      case GLP_FEAS: cout << "solution is feasible;" << endl; break;
      case GLP_INFEAS: cout << "solution is infeasible;" << endl; break;
      case GLP_NOFEAS: cout << "problem has no feasible solution;" << endl; break;
      case GLP_UNBND: cout << "problem has unbounded solution;" << endl; break;
      case GLP_UNDEF: cout << "solution is undefined." << endl; break;
    }

    bool solution_found = glp_prim_stat == GLP_OPT || glp_prim_stat == GLP_FEAS;

    if (solution_found)
    {
      for (int i = 1; i <= variable_size; i++)
      {
        double x = glp_mip_col_val(lp, i);
        if(x!=0)
        {
          //cout<<"solution i="<<i-1<<" index2eid[i]="<<index2eid[i-1]<<endl;
          cuts.push_back( index2eid[i-1] );
        }
      }
    }

    glp_delete_prob(lp);
    delete[] ia;
    delete[] ja;
    delete[] ar;

    return solution_found;
}

//
//
// TopologicalNetSurgent
//
//

TopologicalNetSurgent::TopologicalNetSurgent(int max_iter)
{
  this->m_max_iterations = max_iter;
}

bool compareNetTree(TopologicalNetSurgent::NetTree * t1, TopologicalNetSurgent::NetTree * t2)
{
  return t1->getError()<t2->getError();
}

//convert the net into a net set merge them into the fewest number of
//disjoint nets
NetSet * TopologicalNetSurgent::apply(Net * net)
{
  //get an initial netset by splitting the net
  NetSet * netset = SetCoverNetSurgent::apply(net);

  NetGraph netG(netset);

  vector<NetMatch> all_matches;
  uint size=netG.getVertexSize();
  cout<<"- [Topo Net Surgent] There are "<<size<<" nets to merge"<<endl;
  //uint count=0;
  for(uint i=0;i<size;i++)
  {
    for(uint j=i+1;j<size; j++)
    {
      //cout<<"i="<<i<<" j="<<j<<endl;
#if 1 //use only the best
      NetMatch match=getBestMatch(netset, netG.getVertex(i), netG.getVertex(j));
      if(match.error==DBL_MAX) continue; //invalid matching
      all_matches.push_back(match);
      netG.addEdge(i,j,match.eid);
#else
      vector<NetMatch> tmp;
      getMatches(netset, netG.getVertex(i), netG.getVertex(j), tmp);
      for(NetMatch& match : tmp)
      {
        netG.addEdge(i,j,match.eid);
      }
      all_matches.insert(all_matches.end(),tmp.begin(),tmp.end());
#endif //

      //count++;
    }//end j
  }//end i

  //create some spanning trees and evaluate their nets
  const double original_error = NetSet(net).evaluate();
  const int generation=20;
  int population=size*2;

  cout<<"- [Topo Net Surgent] original_error="<<original_error<<endl;
  vector<TopologicalNetSurgent::NetTree *> trees;
  cout<<"- [Topo Net Surgent] create initial population of "<<population<<" "<<flush;
  int failed_count=0;
  for(int i=0;i<population && failed_count<100 ;i++)
  {
      auto tree=netG.getRandomST();
      if(tree->getError()>=original_error){ i--; failed_count++; continue; }
      trees.push_back(tree);
      //cout<<"\ttree error="<<tree->getError()<<endl;
      cout<<"."<<flush;
      failed_count=0;
  }
  cout<<"(done)"<<endl;
  population=trees.size();
  if(population==0) return netset; //nothing to evolve

  cout<<"- [Topo Net Surgent] start evolving "<<flush;
  for(int i=0;i<generation;i++)
  {
      //create new generations by mutation
      for(int j=0;j<population/2;j++)
      {
        //cout<<"Before mutate"<<endl;
        uint r=j;//(uint)floor(drand48()*population);
        trees.push_back(trees[r]->Mutate());
        //cout<<"after mutation "<<trees.back()->getError()<<endl;
        cout<<"."<<flush;
      }

      //create new generations by crossing over
      for(int j=0;j<population/2;j++)
      {
        uint s=(uint)floor(drand48()*population);
        uint t=s;
        while(t==s) t=(uint)floor(drand48()*population);
        trees.push_back(trees[s]->Crossover(trees[t]));
        //cout<<"after cross over "<<trees.back()->getError()<<endl;
        cout<<"."<<flush;
      }

      //sort by error, smallest error first
      sort(trees.begin(), trees.end(), compareNetTree );

      //remove second half
      //cout<<"remove second half"<<endl;
      //cout<<"trees size="<<trees.size()<<endl;
      for(int j=population;j<trees.size();j++){
        delete trees[j];
      }
      trees.resize(population);

      cout<<"\n- [Topo Net Surgent] [generation "<< i <<"] lowest error="<<trees.front()->getError()<<flush;
  }
  cout<<endl;

  //return the first
  if(trees.front()->getError()<original_error)
  {
    delete netset;
    netset=trees.front()->getNetSet();
    assert(netset->getNets().size()==1);
    netset = apply(netset->getNets().front());
  }
  else
  {
    cout<<"- [Topo Net Surgent] No better netset is found. Done."<<endl;
  }

  //free
  for(int j=0;j<trees.size();j++){
    if(netset!=trees[j]->getNetSet()) delete trees[j];
  }

  return netset;
}

TopologicalNetSurgent::NetMatch
TopologicalNetSurgent::getBestMatch(NetSet * netset, Net * n1, Net * n2)
{
  NetMatch best_match(n1,n2);
  set<uint> sharedE;
  netset->sharedCutEdges(n1, n2, sharedE);
  if(sharedE.empty()) return best_match;

  //cout<<"========================\n";//"\n n1="<<*n1<<"\n n2="<<*n2<<endl;
  for(uint eid : sharedE)
  {
    double error=netset->evaluate_merge(eid);
    //cout<<"error = "<<error<<" by merge nets at eid="<<eid<<endl;
    if(error<best_match.error)
    {
      best_match.error=error;
      best_match.eid=eid;
    }
  }

  return best_match;
}


void
TopologicalNetSurgent::getMatches
(NetSet * netset, Net * n1, Net * n2, vector<TopologicalNetSurgent::NetMatch>& matches)
{
  set<uint> sharedE;
  netset->sharedCutEdges(n1, n2, sharedE);

  //cout<<"========================\n";//"\n n1="<<*n1<<"\n n2="<<*n2<<endl;
  for(uint eid : sharedE)
  {
    NetMatch match(n1,n2);
    match.error=netset->evaluate_merge(eid);
    match.eid=eid;
    if(match.error==DBL_MAX) continue; //invalid matching
    matches.push_back(match);
  }

}

TopologicalNetSurgent::NetGraph::NetGraph(NetSet * netset)
{
  this->m_netset = netset;
  const auto& nets = netset->getNets();
  m_nets=vector<Net*>(nets.begin(),nets.end());
  int vsize=m_nets.size();
  this->resize(vsize);
  for(int i=0;i<vsize;i++) (*this)[i].resize(vsize);
  m_is_vertex_pairs_a_set=false;
}

//create a random spanning tree
TopologicalNetSurgent::NetTree *
TopologicalNetSurgent::NetGraph::getRandomST()
{
  if(!m_is_vertex_pairs_a_set)
  {
    //convert m_is_vertex_pairs_a_set to a set
    set< pair<uint,uint> > tmp(m_vertex_pairs.begin(),m_vertex_pairs.end());
    m_vertex_pairs=vector< pair<uint,uint> >(tmp.begin(),tmp.end());
    m_is_vertex_pairs_a_set=true;
  }

  uint vsize=this->getVertexSize();
  uint pairsize=this->m_vertex_pairs.size();

  DisjointSets ds(vsize);
  NetTree * tree=new NetTree(this);
  random_shuffle(m_vertex_pairs.begin(),m_vertex_pairs.end());

  int count=0;
  for(const auto & vpair : m_vertex_pairs)
  {
    uint s=vpair.first;
    uint t=vpair.second;
    if(ds.find(t)==ds.find(s)) continue; //already connected
    int esize=(*this)[s][t].size();
    if(esize==0){
      cerr<<"! Error: TopologicalNetSurgent::NetGraph::getRandomMST Error"<<endl;
      continue; //no edge?? this should not happend
    }
    uint e=(uint)floor(drand48()*esize);
    if( tree->addEdge(s,t,(*this)[s][t][e])==false ) continue;
    ds.unite( ds.find(t), ds.find(s) );
    count++;
    if(count==(vsize-1)) break; //no more than (vsize-1) edge
  }

  //cout<<"tree netset size="<<tree->getNetSet()->getNets().size()<<" count="<<count<<endl;

  return tree;
}

TopologicalNetSurgent::NetTree::NetTree(NetGraph * source)
{
  this->m_source = source;
  int vsize=source->size();
  this->resize(vsize);
  for(int i=0;i<vsize;i++) (*this)[i].resize(vsize, UINT_MAX);

  //clonse the netset from source...
  this->m_netset=new NetSet(source->m_netset);

  //init
  error=DBL_MAX;
}

TopologicalNetSurgent::NetTree::~NetTree()
{
  delete this->m_netset;
}

bool TopologicalNetSurgent::NetTree::addEdge(uint s, uint t, uint e)
{
  if(this->m_netset->merge(e)==false){
    cerr<<"! Error: TopologicalNetSurgent::NetTree::addEdge Failed"<<endl;
    return false;
  }

  (*this)[s][t]=(*this)[t][s]=e;
  //also merge nets in the netset
  m_vertex_pairs.push_back((s<t)?make_pair(s,t):make_pair(t,s) );
  return true;
}

TopologicalNetSurgent::NetTree * TopologicalNetSurgent::NetTree::Mutate()
{
  //pick 1/2 edges from this and pick the rest randomly
  //cout<<"from tree"<<endl;
  uint vsize=this->m_source->getVertexSize();
  DisjointSets ds(vsize);
  uint esize=m_vertex_pairs.size();
  NetTree * tree=new NetTree(this->m_source);
  uint from_this=(uint)floor(esize*0.5);
  random_shuffle(m_vertex_pairs.begin(),m_vertex_pairs.end());
  for(uint i=0;i<from_this;i++){
    const auto & vpair = m_vertex_pairs[i];
    uint s=vpair.first;
    uint t=vpair.second;
    if(tree->addEdge(s,t,(*this)[s][t])==false){ i--; continue; }
    ds.unite( ds.find(t), ds.find(s) );
  }

  //cout<<"from graph"<<endl;
  //add the rest randomly using edges in this->m_source
  const auto& source_vpairs=this->m_source->m_vertex_pairs;
  int totalesize=source_vpairs.size();
  for(uint i=from_this;i<esize;i++)
  {
    //get a random pair from source_vpairs
    const auto & vpair = source_vpairs[(uint)floor(drand48()*totalesize)];
    uint s=vpair.first;
    uint t=vpair.second;
    if(ds.find(t)==ds.find(s)){ i--; continue; }//already connected
    int esize=(*m_source)[s][t].size();
    if(esize==0){
      cerr<<"! Error: TopologicalNetSurgent::NetTree::Mutate Error"<<endl;
      continue; //no edge?? this should not happend
    }
    uint e=(uint)floor(drand48()*esize);
    if( tree->addEdge(s,t,(*m_source)[s][t][e])==false ){ i--; continue; }
    ds.unite( ds.find(t), ds.find(s) );
  }//end for i

  return tree;
}

TopologicalNetSurgent::NetTree * TopologicalNetSurgent::NetTree::Crossover(NetTree * other)
{
  //create a new tree
  NetTree * tree=new NetTree(this->m_source);

  //get all edge and mix together
  set< pair<uint,uint> > tmp(m_vertex_pairs.begin(),m_vertex_pairs.end());
  tmp.insert(other->m_vertex_pairs.begin(),other->m_vertex_pairs.end());
  vector< pair<uint,uint> > all_vpairs(tmp.begin(),tmp.end());
  random_shuffle(all_vpairs.begin(),all_vpairs.end());

  //get disjoint set
  uint vsize=this->m_source->getVertexSize();
  DisjointSets ds(vsize);

  //start adding edges
  uint esize=m_vertex_pairs.size();
  int count=0;
  for(const auto & vpair : all_vpairs)
  {
    uint s=vpair.first;
    uint t=vpair.second;
    if(ds.find(t)==ds.find(s)) continue; //already connected
    //cout<<"cross over add edge["<<s<<"]["<<t<<"]="<<(*this)[s][t]<<endl;
    uint e1=(*this)[s][t];
    uint e2=(*other)[s][t];
    uint e=0;
    if(e1!=UINT_MAX && e2!=UINT_MAX)
    {
      e=(drand48()>0.5)?e1:e2;
    }
    else if(e1!=UINT_MAX){ e=e1; }
    else if(e2!=UINT_MAX){ e=e2; }
    else {
      cerr<<"! Error: TopologicalNetSurgent::NetTree::Crossover"<<endl;
      continue; //something is wrong...
    }
    if( tree->addEdge(s,t,e)==false ) continue;
    ds.unite( ds.find(t), ds.find(s) );
    count++;
    if(count==esize) break; //no more than esize edge
  }

  return tree;
}

double TopologicalNetSurgent::NetTree::getError()
{
  if(this->error==DBL_MAX)
  {
    this->error=this->m_netset->evaluate();
  }

  return this->error;
}

//
BoxingNetSurgent::BoxingNetSurgent(float width, float height)
{
  m_box_width=width;
  m_box_height=height;
}


//convert the net into a net set
NetSet * BoxingNetSurgent::apply(Net * net)
{

  cout<<"- Apply BoxingNetSurgent"<<endl;

  if(valid(net)) //check if the input net is valid
  {
    return new NetSet(net);
  }

  //NetSet * netset = new NetSet(net);
  polygon::contained_bbox problem(m_box_width, m_box_height);

  vector<_NetSetCandidate> open;
  open.push_back(_NetSetCandidate(net));
  NetSet* solution=NULL;
  list<uint> solution_cut_edges;
  bool valid_net_found=false;

  while(open.empty()==false)
  {
    _NetSetCandidate ns=open.front();
    pop_heap(open.begin(),open.end());
    open.pop_back();

    //check if we find a solution
    if(ns.invalid_net==NULL || ns.invalid_net->getFaces().size()<3)
    {
      //solution=ns.ns;
      solution_cut_edges=ns.cut_edges;
      break;
    }

    {
      auto poly=ns.invalid_net->getNetBoundary();
    //   cout<<"polygon=\n"<<poly<<endl;
      polygon::min_perimeter_bbox problem;
      polygon::bbox2d solver(poly);
      polygon::obb box = solver.build(problem);
      cout<<"bbox="<<box<<endl;
    }

    //
    cout<<"- Open size = "<<open.size()<<", Solution "
        <<((valid_net_found)?"found    ":"not found")
        <<"\r"<<flush;
    //cout<<"===GET 10 BEST============== open size="<<open.size()<<endl;
    vector<_NetSetCandidate> open2;
    //cout<<"netset size="<<ns.ns->getNets().size()<<endl;
    for(uint eid : ns.invalid_net->getCreases())
    {
      //cout<<"-----------------------------"<<endl;
      //cout<<"eid="<<eid<<endl;
      pair<Net *, Net *> new_nets=ns.invalid_net->split(eid);

      //sometime wrong, ignore the edge
      if(new_nets.second==NULL) continue;

      //nets are too small to be interesting
      if(new_nets.first->getFaces().size()<ns.invalid_net->getFaces().size()/10 ||
         new_nets.second->getFaces().size()<ns.invalid_net->getFaces().size()/10)
      {
        if(new_nets.first!=NULL) delete new_nets.first;
        if(new_nets.second!=NULL) delete new_nets.second;
        continue;
      }

      //cout<<"net 1 has "<<new_nets.first->getFaces().size()<<" faces"<<endl;
      bool r1=valid(new_nets.first);
      //cout<<"net 2 has "<<new_nets.second->getFaces().size()<<" faces"<<endl;
      bool r2=valid(new_nets.second);

      if( r1||r2 ) //one of the nets is valid
      {
        Net * invalid_net = (!r1)?new_nets.first:((!r2?new_nets.second:NULL));
        _NetSetCandidate _nsc(invalid_net);
        _nsc.eid=eid;
        _nsc.cut_edges=ns.cut_edges;
        _nsc.cut_edges.push_back(eid);
        _nsc.net_size=ns.net_size+1;
        open2.push_back(_nsc);
        push_heap(open2.begin(),open2.end());
        if(_nsc.invalid_net==NULL) valid_net_found=true;
      }
      else //both too big
      {
        //delete new_nets.first;
        //delete new_nets.second;
      }
    }//end for eid

    for(int i=0;i<2;i++) //only keep 10 best and insert them to open
    {
      open.push_back(open2.front());
      //if(open2.front().invalid_net!=NULL)
      //cout<<"open2.front() area="<<open2.front().max_bbox.width*open2.front().max_bbox.height<<endl;
      //cout<<"eid="<<open2.front().eid<<endl;
      push_heap(open.begin(),open.end());
      pop_heap(open2.begin(),open2.end());
      open2.pop_back();
      if(open2.empty()) break;
    }

    //done with this ns
    //delete ns.ns;

    //break;

  }//end while

  //release memory
  // for(_NetSetCandidate& ns : open)
  // {
  //     delete ns.ns;
  // }

  // cout<<"cut edges=";
  // for(uint eid : solution_cut_edges) cout<<eid<<" ";
  // cout<<endl;

  if(!solution_cut_edges.empty())
  {
    solution = new NetSet(net);
    //cout<<"SPLIT 1="<<solution->split(solution_cut_edges.back())<<endl;
    //cout<<"SPLIT 2="<<solution->split(solution_cut_edges.front())<<endl;;
    for(uint eid : solution_cut_edges){
      solution->split(eid);
      //break;
    }
  }

  // cout<<"orig net size="<<net->getFaces().size()<<endl;
  // for(Net * net : solution->getNets())
  // {
  //     cout<<"solution net\n="<<net->getNetBoundary().front()<<endl;
  // }
  cout<<endl;
  return solution;
}

BoxingNetSurgent::_NetSetCandidate::_NetSetCandidate(Net * _net)
{
  invalid_net=_net;
  net_size=1;
  if(invalid_net!=NULL)
  {
    const polygon::c_polygon & p=invalid_net->getNetBoundary();
    polygon::min_perimeter_bbox problem;
    polygon::bbox2d solver(p);
    this->max_bbox = solver.build(problem);
  }
}

//compare two netsets
bool BoxingNetSurgent::_NetSetCandidate::operator<(const _NetSetCandidate& other) const
{
  //cout<<"WHATS???"<<endl;
  if(this->invalid_net==NULL && other.invalid_net!=NULL) //this is a valid net
  {
    //cout<<"YESSSS"<<endl;
    return false;
  }

  if(this->invalid_net!=NULL && other.invalid_net==NULL) //other is a valid net
  {
    //cout<<"NOOOOOO"<<endl;
    return true;
  }

  //both nets are invalid
  if(this->net_size==other.net_size)
  {
    //cout<<"YEAH!!!!!"<<endl;

    //else (both null or both non-null)
    auto& b1=this->max_bbox;
    auto& b2=other.max_bbox;

    return (b1.width+b1.height)>(b2.width+b2.height);
  }

  //cout<<"!!!!!!!!"<<endl;
  return this->net_size>other.net_size;
}

//
bool BoxingNetSurgent::valid(Net * net)
{
  if(net==NULL) return false; //safe guard

  polygon::contained_bbox problem(m_box_width, m_box_height);

  if( m_valid_nets.find(net) == m_valid_nets.end() )
  {
    //compute the bounding box of the net and check if it is valid
    //auto poly=net->getNetBoundary();
    //cout<<"net polygon="<<poly<<endl;

    polygon::bbox2d solver(net->getNetBoundary());
    polygon::obb box = solver.build(problem);

    //cout<<box<<endl;

    //if not valid, return FALSE
    //else save the net into the set
    if(problem.solved(box))
    {
      m_valid_nets.insert(net);
    }
    else return false;
  }

  return true;
}

//
bool BoxingNetSurgent::valid(NetSet * netset)
{
  for(Net * net : netset->getNets())
    if(valid(net)==false) return false;
  return true;
}

//each netset must only have one invalid net
Net * BoxingNetSurgent::getInvalidNet(NetSet * netset)
{
  for(Net * net : netset->getNets())
    if(valid(net)==false) return net;
  return NULL;
}

}//end namespace masc
