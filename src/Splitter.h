/*
 * Splitter.h
 *
 * A list of classes that computed edge weights for unfolding
 * Will be better if this is called a "Weighter" as these classes does not split
 * anything....
 *
 *  Created on: Feb 6, 2015
 *      Author: zhonghua
 */

#ifndef SPLITTER_H_
#define SPLITTER_H_

#include <vector>
using namespace std;

#include "config.h"

struct model;
class Unfolder;

namespace masc {

class Splitter {
public:
    virtual void measure(model* m);
    virtual ~Splitter(){};
    // Assign dual edge weights, weights should be in [0, 1]
    vector<float> assignWeights(model* m, const Config& config);
    static Splitter* createSplitter(CutHeuristic heruistic);
    static Splitter* createSplitter(const string& heruistic);

    virtual string getName()=0;

protected:
    void virtual assignWeightsImpl(model* m, vector<float>& weights, const Config& config)=0;
    // normalize the given weights.
    void normalizeWeights(vector<float>& weights);
    // generate a random unit vector
    virtual  Vector3d genRandomUnitVector(const Config& config);
    Splitter();
    float m_min_edge_length;
    float m_max_edge_length;
    float m_total_edge_length;
private:
    void trimWeights(vector<float>& weights);
};

class FlatTreeSplitter : public Splitter {
public:
    FlatTreeSplitter(){};
    virtual ~FlatTreeSplitter(){};
    string getName() override { return "FlatTreeSplitter"; }
protected:
    void assignWeightsImpl(model* m, vector<float>& weights, const Config& config) override;
};

class UnFlatTreeSplitter : public FlatTreeSplitter {
public:
    UnFlatTreeSplitter(){};
    virtual ~UnFlatTreeSplitter(){};
    string getName() override { return "UnFlatTreeSplitter"; }
protected:
    void assignWeightsImpl(model* m, vector<float>& weights, const Config& config) override;
};

// to minimize the perimeter of the unfolding
class MinimumPerimeterSplitter : public Splitter {
public:
    MinimumPerimeterSplitter(){}
    virtual ~MinimumPerimeterSplitter(){}
    string getName() override { return "MinimumPerimeterSplitter"; }
protected:
    void assignWeightsImpl(model* m, vector<float>& weights, const Config& config) override;
};

class MaximumPerimeterSplitter : public MinimumPerimeterSplitter {
public:
    MaximumPerimeterSplitter(){}
    virtual ~MaximumPerimeterSplitter(){}
    string getName() override { return "MaximumPerimeterSplitter"; }
protected:
    void assignWeightsImpl(model* m, vector<float>& weights, const Config& config) override;
};

class SteepestEdgeSplitter : public Splitter {
public:
    SteepestEdgeSplitter() : m_steepest(true) {}
    virtual ~SteepestEdgeSplitter() {}
    string getName() override { return "SteepestEdgeSplitter"; }
protected:
    void assignWeightsImpl(model* m, vector<float>& weights, const Config& config) override;
    bool m_steepest;
};

class UnSteepestEdgeSplitter : public SteepestEdgeSplitter {
public:
    UnSteepestEdgeSplitter() { m_steepest = false; }
    virtual ~UnSteepestEdgeSplitter() {}
    string getName() override { return "UnSteepestEdgeSplitter"; }
};

// assign random weights on edges
class RandomSplitter : public Splitter {
public:
    RandomSplitter(){}
    virtual ~RandomSplitter() {}
    string getName() override { return "RandomSplitter"; }
protected:
    void assignWeightsImpl(model* m, vector<float>& weights, const Config& config) override;
};

// assign random weights on edges for each vertex
class RandomEdgeSplitter : public Splitter {
public:
    RandomEdgeSplitter(){}
    virtual ~RandomEdgeSplitter() {}
    string getName() override { return "RandomEdgeSplitter"; }
protected:
    void assignWeightsImpl(model* m, vector<float>& weights, const Config& config) override;
};

// try all possible random tree
class BruteForceSplitter : public Splitter {
public:
    BruteForceSplitter():m_inited(false) {}
    virtual ~BruteForceSplitter() {}
    string getName() override { return "BruteForceSplitter"; }
protected:
    void assignWeightsImpl(model* m, vector<float>& weights, const Config& config) override;
    void init(int edges);

    vector<float> m_weights;
    bool m_inited;
};

} /* namespace masc */

#endif /* SPLITTER_H_ */
