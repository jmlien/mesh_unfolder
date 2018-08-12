#include "netsurgent.h"
#include "unfolder.h"
#include <unordered_map>
#include <deque>

namespace masc {

bool NetSurgent::operate(Unfolder* unfolder, vector<Unfolder*> & operated)
{
  if (unfolder->isFlattened()) return false; //no need to operate

  //get config
  const Config& config = unfolder->getConfig();

  //get model and creases to create a net
  cout<<"create net"<<endl;
  model * m=unfolder->getModel();
  const set<uint>& creases=unfolder->getFoldEdges();
  Net * net=new Net(m,creases);

  //create a net surgent and apply
  cout<<"create surgent"<<endl;
  NetSurgent * surgent=NULL;
  switch(config.surgery_method)
  {
    case NetSurgery::SET_COVER_SURGERY: surgent=new SetCoverNetSurgent(); break;
    case NetSurgery::TOPOLOGICAL_SURGERY:
    case NetSurgery::SUBDIVID_SURGERY:
    case NetSurgery::CAGING_SURGERY:
    default:
      cerr<<"! Error: Net Surgery type unsupported"<<endl;
  }

  if(surgent==NULL) return false;

  cout<<"surgent operates the net"<<endl;
  NetSet * netset = surgent->apply(net);
  list<Net*> & nets=netset->getNets();

  //convert nets to a list of new unfolders
  cout<<"convert nets to unfolders"<<endl;
  int cid=1;
  for(Net* net: nets)
  {
    Unfolder * tmp=net->toUnfolder(config);
    tmp->setClusterId(cid++);
    operated.push_back(tmp);
  }

  //done
  return true;
}

SetCoverNetSurgent::SetCoverNetSurgent()
{

}

//convert the net into a net set
NetSet * SetCoverNetSurgent::apply(Net * net)
{
    //find all all_overlaps
    set< pair<uint,uint> >  overlaps = net->getOverlaps();

    cout<<"there are "<<overlaps.size()<<" overlaps"<<endl;

    //find shortest paths (sequence of edges) connecting pairs of overlapping Triangles
    unordered_map<uint, list<pair<uint,uint> > > e2tmap; //map< eid, pair<uint, uint>=<tid1, tid2> >
    for(const pair<uint,uint>& overlap : overlaps)
    {
      list<uint> epath;
      cout<<"find path for "<<overlap.first<<","<<overlap.second<<": ";
      //find path between overlap.first and overlap.second
      this->find_path(net, overlap.first, overlap.second, epath);

      for(uint e : epath)
      {
        e2tmap[e].push_back(overlap);
        cout<<e<<", ";
      }
      cout<<endl;
    }

    cout<<"finished finding all paths, e2tmap has size="<<e2tmap.size()<<endl;

    //this might be slow
    //for each edge, we maintain a list of overlapping pair
    //greedily, we start from an edge with the largest number
    //of un-resolved pairs
    cout<<"solve set-cover problem"<<endl;
    list<uint> new_cuts;
    while(overlaps.empty()==false)
    {
      //find an edge with the largest cover
      uint best_e;
      uint best_size=0;
      for(auto& tmp : e2tmap) //loop through each edge
      {
          uint eid = tmp.first;
          const list<pair<uint,uint> >& conflicts = tmp.second;
          //count the number of resolved pairs
          int count=0;
          for(const pair<uint,uint> & overlap : conflicts)
          {
            if(overlaps.find(overlap)!=overlaps.end()) count++;
          }
          //
          if(count>best_size)
          {
            best_size=count;
            best_e=eid;
          }
          else if(count==best_size) //resolve the same number of pairs
          {
            //pick one with shorter length
            model * m = net->getOriginalModel();
            edge & e1 = m->edges[best_e];
            edge & e2 = m->edges[eid];
            if(e2.length<e1.length)
            {
              best_e=eid;
            }
          }
      }//end for tmp

      cout<<"best_e="<<best_e<<" best_size="<<best_size<<endl;

      //found the best edge, remember it and remove all resolved pairs from overlaps
      cout<<"remove resolved overlaps"<<endl;
      new_cuts.push_back(best_e);
      for(const pair<uint,uint> & overlap : e2tmap[best_e])
      {
        overlaps.erase(overlap); //these overlaps have been resolved
      }//end for overlap


      cout<<"remove best_e"<<endl;
      e2tmap.erase(best_e); //remove best_e from the map

      cout<<"overlaps size="<<overlaps.size()<<endl;
    }//end while

    //cut the edges
    cout<<"cut nets"<<endl;
    NetSet * netset = new NetSet(net);
    for(uint eid : new_cuts)
    {
      cout<<"cut edge "<<eid<<endl;
      bool r = netset->split(eid);
      if(!r) cerr<<"! Error: failed to cut edge "<<eid<<endl;
      //break;
    }

    cout<<"surgery done"<<endl;
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


}//end namespace masc
