#include "net.h"

namespace masc {

//
//a set of nets spanning the given mesh
//
class NetSet{

public:

    NetSet(Net * n);

    virtual ~NetSet();

    bool split(uint eid); //split a net at a given edge, this edge must be a crease edge

    bool merge(uint eid); //merge two nets at a given edge, this edge must be a cut edge

    list<Net*> & getNets(){ return m_nets; }

protected:

  void replace(Net *net);

private:

    model * m_orig;

    list<Net*> m_nets;

    vector<Net*> m_fid2net; //reverse query, given an fid, record the net that contains this face
};

} //namespace masc
