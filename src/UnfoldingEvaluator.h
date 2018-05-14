/*
 * UnfoldingEvaluator.h
 *
 *  Created on: Feb 29, 2016
 *      Author: zxi
 */

#ifndef SRC_UNFOLDINGEVALUATOR_H_
#define SRC_UNFOLDINGEVALUATOR_H_

#include <memory>
#include <map>
#include <string>

namespace masc {
namespace ga {
// forward declaration
class Individual;
}
}

class Unfolder;

namespace masc {
namespace unfolding {

// forward declaration
class OverlappingChecker;

///////////////////////////////////////////////////////////////////////////////
// Base Evaluators
///////////////////////////////////////////////////////////////////////////////

// Evaluator: measure the goodness of an unfolding
class UnfoldingEvaluator {
public:
  // evaluate an individual
  UnfoldingEvaluator(Unfolder* unfolder) :
      m_unfolder(unfolder) {
  }
  virtual ~UnfoldingEvaluator() {
    m_unfolder = nullptr;
  }
  virtual double evaluate(masc::ga::Individual* ind) = 0;
protected:
  // does not own the object
  Unfolder* m_unfolder;
};

// Evaluate the fitness of the net after flatten
class NetEvaluator {
public:
  virtual ~NetEvaluator() {}

  // Should return a number > 0, the larger the better
  virtual double evaluate(Unfolder* unfolder) = 0;
};

///////////////////////////////////////////////////////////////////////////////
// OverlappingEvaluator
//
// Measure the goodness as number of overlapping
///////////////////////////////////////////////////////////////////////////////

class OverlappingEvaluator: public UnfoldingEvaluator {
public:
  OverlappingEvaluator(Unfolder* unfolder) :
      UnfoldingEvaluator(unfolder) {
  }
  virtual ~OverlappingEvaluator() {
  }
  virtual double evaluate(masc::ga::Individual* ind) override;
};



class CutLengthEvaluator : public NetEvaluator {
public:
  virtual ~CutLengthEvaluator(){}
  double evaluate(Unfolder* unfolder) override;
};

class HullAreaEvaluator : public NetEvaluator {
public:
  virtual ~HullAreaEvaluator(){}
  double evaluate(Unfolder* unfolder) override;
};

// estimate max side length of bounding box
class BoxSideEvaluator : public NetEvaluator {
public:
  virtual ~BoxSideEvaluator(){}
  double evaluate(Unfolder* unfolder) override;
};

//evaluate based how similar the net bounrady is to a list of pre-define polygons
class PolygonFitEvaluator : public NetEvaluator {
public:

  PolygonFitEvaluator(const std::string& stencil_filename);
  virtual ~PolygonFitEvaluator();
  double evaluate(Unfolder* unfolder) override;

private:
  void * m_poly_ptr;
  void * m_best_net_ptr;
  double m_min_error;
  int best_src_off, best_src_len, best_target_off, best_target_len;
};

///////////////////////////////////////////////////////////////////////////////
// LearningEvaluator
//
// Learn a fitness function
///////////////////////////////////////////////////////////////////////////////

class LearningEvaluator: public UnfoldingEvaluator {
public:
LearningEvaluator(Unfolder* unfolder, std::map<std::string,std::string> * params) :
      UnfoldingEvaluator(unfolder) {
          m_params = params;
  }
  virtual ~LearningEvaluator() {
  }
  virtual double evaluate(masc::ga::Individual* ind) override;
protected:
  std::map<std::string, std::string> *m_params;
};

///////////////////////////////////////////////////////////////////////////////
// AreaEvaluator
//
// Measure the goodness as area ratio
///////////////////////////////////////////////////////////////////////////////

class AreaEvaluator: public UnfoldingEvaluator {
public:
  AreaEvaluator(Unfolder* unfolder);
  virtual ~AreaEvaluator() {
  }
  virtual double evaluate(masc::ga::Individual* ind) override;
protected:
  std::unique_ptr<OverlappingChecker> m_checker;
  double m_best_ratio;
};

} /* namespace unfoldings */
} /* namespace masc */

#endif /* SRC_UNFOLDINGEVALUATOR_H_ */
