//this class represents a net structure that can be merged with other nets, split into multiple nets
//
//

#include "model.h"
#include "polygon/bbox2d.h"

class Unfolder; //defined in unfolder.h
struct Config;  //defined in config.h

namespace masc {

class Net
{
public:

    Net(model * m, const set<uint>& crease_edges, Net * p1=NULL, Net * p2=NULL);
    Net(Net * n); //clone a net...

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
    const set<uint>& getCreases() const {return m_creases; }

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

    const polygon::c_polygon & getNetBoundary() const { return m_boundary; }

protected:

    //create a net with a single face...
    Net(model * m, int fid, Net * p1);

    //build subm from m_orig and m_faces
    void build_subm();

    //unfold m_subm into m_unfolding
    void unfold_subm();

private:

    //split the net using a seed_fid and cut edge id (eid)
    Net * split(uint eid, uint seed_fid, vector<bool>& face_visited);

    //data
    Net * m_parents[2]; //where this net is from

    //these are local data and ids of m_subm
    shared_ptr<model> m_subm; //subm for this net

    //MESH  m_subm_unfolded; //unfolding of the m_subm
    vector<bool> m_subm_crease_lines; //in subm, which edges are crease lines.

    //these are data and ids of m_orig
    model * m_orig; //original model that this net comes from
    list<uint> m_faces;
    set<uint> m_creases;
    set< pair<uint,uint> > m_overlap_pairs; //pairs of overlapping faces

    //boundary of this net
    polygon::c_polygon m_boundary;
};


//
inline ostream & operator<<(ostream & out, Net & net) {
  out<<"(";
  for(auto& id : net.getFaces()) out<<id<<", ";
  out<<")";
  return out;
}


}//end namespace masc
