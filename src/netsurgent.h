#include "netset.h"
#include <cfloat>
#include <climits>

class Unfolder; //defined in unfolder.h

namespace masc {


class NetSurgent
{
public:

    static bool operate(Unfolder* unfolder, vector<Unfolder*>& operated);

    //convert the net into a net set
    virtual NetSet * apply(Net * net)=0;
};

//this is class that converts a net with overlaps into a
//set of non-overlapping nets using set cover
class SetCoverNetSurgent : public NetSurgent
{
public:

    SetCoverNetSurgent();

    //convert the net into a net set
    virtual NetSet * apply(Net * net);

protected:

    NetSet * apply2(Net * net);

    struct LP_constraints
    {
      LP_constraints(){ type = 1;  upper_bound = lower_bound = 0; }
      vector<uint> eids; //edges involed in this constraint

      int type;  //GLP_FR    free (unbounded) variable, (1)
                 //GLP_LO    lower bound
                 //GLP_UP    upper bound
                 //GLP_DB    double bound
                 //GLP_FX    fixed

      int upper_bound;
      int lower_bound;
    };

    //find the edge-path between f1 and f2 in the net
    void find_path(Net * net, uint f1, uint f2, list<uint>& path);
    bool SolveSetCover(const list< list<uint> > & epaths, list<uint>& cuts);
};//end SetCoverNetSurgent

//
// using merge and split operators to reduce the number of non-overlapping nets
//
class TopologicalNetSurgent : public SetCoverNetSurgent
{
public:

  TopologicalNetSurgent(int max_iter);

  //convert the net into a net set
  virtual NetSet * apply(Net * net);

private:

  struct NetMatch
  {
    NetMatch(Net * n1, Net * n2)
    {
      net1=n1; net2=n2; eid=UINT_MAX; error=DBL_MAX;
    }
    Net * net1, * net2;
    uint eid;
    double error;
  };

  //
  // the idea here is to connect all nets into a graph, whose nodes are Nets
  // and edges are NetMatch
  // Then, we extract many spanning trees to evolve then into a spanning tree with
  // minimum errors
  class NetTree; //defined below

  class NetGraph : public vector< vector< vector<uint> > >
  {
    public:

      NetGraph(NetSet * m_netset);
      NetTree * getRandomST(); //create a random spanning tree
      void addEdge(uint s, uint t, uint e){
        (*this)[s][t].push_back(e);
        (*this)[t][s].push_back(e);
        m_vertex_pairs.push_back( (s<t)?make_pair(s,t):make_pair(t,s) );
      }

      uint getVertexSize(){ return m_nets.size(); }

      Net * getVertex(int i){ return m_nets[i]; }

    private:

      friend class NetTree;
      NetSet * m_netset;
      vector<Net*> m_nets;
      vector< pair<uint,uint> > m_vertex_pairs; //connected vertex pairs
      bool m_is_vertex_pairs_a_set;
  };

  class NetTree : public vector< vector< uint > >
  {
    public:
      NetTree(NetGraph * source);
      virtual ~NetTree();
      NetTree * Mutate();
      NetTree * Crossover(NetTree * other);
      bool addEdge(uint s, uint t, uint e); //return false if failed to merge two nets
      double getError();
      NetSet * getNetSet(){ return m_netset; }

    private:

      NetGraph * m_source;
      NetSet * m_netset;
      double error;
      vector< pair<uint,uint> > m_vertex_pairs; //connected vertex pairs
  };

  friend bool compareNetTree(NetTree * t1, NetTree * t2);

  vector<NetTree*> m_net_pool; //a pool of possible nets

  NetMatch getBestMatch(NetSet * netset, Net * n1, Net * n2);

  void getMatches(NetSet * netset, Net * n1, Net * n2, vector<NetMatch>& matches);

  uint m_max_iterations;//max iterations allowed to improve the net

};//end TopologicalNetSurgent



//
// using merge and split operators to reduce the number of non-overlapping nets
//
class BoxingNetSurgent : public NetSurgent
{
public:

  BoxingNetSurgent(float width, float height);

  //convert the net into a net set
  virtual NetSet * apply(Net * net);

private:

  struct NetAABB{
      Net * net;
      Vector2d u, v; //directions
      float w,h; //size
  };

  NetAABB computeAABB(Net * net);

  float m_box_width, m_box_height; //width and height of the box

};//end BoxingNetSurgent

}//end namespace masc
