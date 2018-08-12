#include "netset.h"

class Unfolder; //defined in unfolder.h

namespace masc {

//this is class that converts a net with overlaps into a
//set of non-overlapping nets using set cover

class NetSurgent
{
public:

    static bool operate(Unfolder* unfolder, vector<Unfolder*>& operated);

    //convert the net into a net set
    virtual NetSet * apply(Net * net)=0;
};

class SetCoverNetSurgent : public NetSurgent
{
public:

    SetCoverNetSurgent();

    //convert the net into a net set
    virtual NetSet * apply(Net * net);

protected:

    //find the edge-path between f1 and f2 in the net
    void find_path(Net * net, uint f1, uint f2, list<uint>& path);

};

}//end namespace masc
