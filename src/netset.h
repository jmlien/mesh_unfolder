#include "net.h"

namespace masc {

//
//a set of nets spanning the given mesh
//
class NetSet{

public:

    NetSet(Net * n);
    NetSet(NetSet * ns);

    virtual ~NetSet();

    bool split(uint eid); //split a net at a given edge, this edge must be a crease edge

    bool merge(uint eid); //merge two nets at a given edge, this edge must be a cut edge

    double evaluate(); //evalute this netset, 0 is the best possible value

    double evaluate_merge(uint eid); //evalute merging two nets at a given edge
                                     //larger return value mean worse merging result
                                     //0 is the best possible value

    list<Net*> & getNets(){ return m_nets; }

    void sharedCutEdges(Net * n1, Net * n2, set<uint>& shared);


protected:

  void replace(Net *net);

private:

    model * m_orig;

    list<Net*> m_nets;

    list<Net*> m_created_nets; //all nets created in this class, for garbage collection

    vector<Net*> m_fid2net; //reverse query, given an fid, record the net that contains this face
};

} //namespace masc
