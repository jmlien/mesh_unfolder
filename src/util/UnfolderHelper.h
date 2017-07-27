/*
 * UnfolderHelper.h
 *
 *  Created on: Mar 4, 2016
 *      Author: zxi
 */

#ifndef SRC_UTIL_UNFOLDERHELPER_H_
#define SRC_UTIL_UNFOLDERHELPER_H_

#include <vector>
#include <unordered_map>

class Unfolder;

namespace masc {
namespace unfolding {
namespace util {

    typedef std::unordered_map<int, std::unordered_map<int, float>> GEO_DIST_GRAPH;
    typedef std::unordered_map<int, std::unordered_map<int, int>> DIST_GRAPH;

class UnfolderHelper {
public:
    UnfolderHelper(Unfolder* unfolder);

  virtual ~UnfolderHelper()
  {

  }

  // compute pairwise face geodesic distances in the unfolding
  GEO_DIST_GRAPH computeGeoDist();

  // compute pairwise face distance in the
  //DIST_GRAPH computeDist();

  // compute diameter of the tree
  int computeDiameter();

  // compute the number of leaf nodes
  int computeNumLeafNodes();

  // compute branch length in the tree
  std::vector<int> computeBranchLength();

  // compute geodesic distances of each branch
  std::vector<float> computeBranchGeoLength();

  // compute face geodesic distances for a given connected component from the root
  std::vector<float> computeGeoDist(const int root, const std::vector<unsigned int>& cc);
  void dumpStats(std::ostream & out);

  void dumpStats(const std::string & path);

  float getTotalCutLength();
  float getTotalCutLengthNormalized();

    std::vector<int> m_num_children;
    float computeVertexDist();
    float computeVertexDistFromSameParent();

    float computeBorderCutsLength();

protected:


  std::vector<int> computeNumChildren();

  // compute geo dist of two
  float computeGeoDist(const int fid1, const int fid2) const;

  // compute geo dist of two
  int computeDist(const int fid1, const int fid2) const;

  // DO NOT own this  object
  Unfolder* m_unfolder;

};

} /* namespace util */
} /* namespace unfolding */
} /* namespace masc */

#endif /* SRC_UTIL_UNFOLDERHELPER_H_ */
