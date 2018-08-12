//this class represents a net structure that can be merged with other nets, split into multiple nets
//
//

#include "model.h"

class Unfolder; //defined in unfolder.h
struct Config;  //defined in config.h

namespace masc {

class Net
{
public:

    Net(model * m, const set<uint>& crease_edges);
    Net(model * m, int fid); //create a net with a single face...

    virtual ~Net();

    //split the net along the edge eid
    pair<Net *, Net *> split(uint eid);

    //merge two nets at a given eid
    Net * merge(Net * n2, uint eid);

    //convert this net into a unfolder
    Unfolder * toUnfolder(const Config & config, bool checkoverlap=false);

    //access
    const list<uint>& getFaces() const { return m_faces; }
    model * getOriginalModel(){ return m_orig; } //return the original model
    const set< pair<uint,uint> > & getOverlaps() { return m_overlap_pairs; }

    //check if a given eid is a crease edge
    //eid is original id in m_orig
    bool is_crease(uint eid)
    {
      //this can be slow...might need to find a constant time solution...
      return m_creases.find(eid)!=m_creases.end();
      // if(eid>m_subm->e_size)
      // {
      //   cerr<<"! Error: is_crease eid out of bound (eid="<<eid<<")"<<endl;
      //   return false;
      // }
      // return m_subm_crease_lines[eid];
    }

protected:

    //build subm from m_orig and m_faces
    void build_subm();

    //unfold m_subm into m_unfolding
    void unfold_subm();

private:

    //split the net using a seed_fid and cut edge id (eid)
    Net * split(uint eid, uint seed_fid, vector<bool>& face_visited);

    model * m_orig; //original model that this net comes from
    model * m_subm; //subm for this net
    MESH  m_unfolded; //unfolding of the m_subm
    set< pair<uint,uint> > m_overlap_pairs; //pairs of overlapping faces
    vector<bool> m_subm_crease_lines; //in subm, which edges are crease lines.

    list<uint> m_faces;
    set<uint> m_creases;
};

}//end namespace masc
