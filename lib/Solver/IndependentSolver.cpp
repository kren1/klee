//===-- IndependentSolver.cpp ---------------------------------------------===//
//
//                     The KLEE Symbolic Virtual Machine
//
// This file is distributed under the University of Illinois Open Source
// License. See LICENSE.TXT for details.
//
//===----------------------------------------------------------------------===//

#define DEBUG_TYPE "independent-solver"
#include "klee/Solver.h"

#include "klee/Expr.h"
#include "klee/Constraints.h"
#include "klee/SolverImpl.h"
#include "klee/Internal/Support/Debug.h"

#include "klee/util/ExprUtil.h"
#include "klee/util/Assignment.h"

#include "llvm/ADT/BitVector.h"
#include "llvm/Support/raw_ostream.h"
#include <list>
#include <map>
#include <ostream>
#include <vector>

using namespace klee;
using namespace llvm;

template <class T>
inline llvm::raw_ostream &operator<<(llvm::raw_ostream &os,
                                     const ::DenseSet<T> &dis) {
  dis.print(os);
  return os;
}

class IndependentElementSet {
public:
  typedef std::map<const Array *, llvm::BitVector> elements_ty;
  elements_ty elements;                 // Represents individual elements of array accesses (arr[1])
  std::set<const Array*> wholeObjects;  // Represents symbolically accessed arrays (arr[x])
  std::vector<ref<Expr> > exprs;        // All expressions that are associated with this factor
                                        // Although order doesn't matter, we use a vector to match
                                        // the ConstraintManager constructor that will eventually
                                        // be invoked.

  IndependentElementSet() = delete;
  IndependentElementSet(ref<Expr> e) {
    exprs.push_back(e);
    // Track all reads in the program.  Determines whether reads are
    // concrete or symbolic.  If they are symbolic, "collapses" array
    // by adding it to wholeObjects.  Otherwise, creates a mapping of
    // the form Map<array, set<index>> which tracks which parts of the
    // array are being accessed.
    std::vector< ref<ReadExpr> > reads;
    findReads(e, /* visitUpdates= */ true, reads);
    for (unsigned i = 0; i != reads.size(); ++i) {
      const ReadExpr *re = reads[i].get();
      const Array *array = re->updates.root;
      
      // Reads of a constant array don't alias.
      if (re->updates.root->isConstantArray() && re->updates.head.isNull())
        continue;

      // Skip rest if already marked as symbolic
      if (wholeObjects.count(array))
        continue;

      if (ConstantExpr *CE = dyn_cast<ConstantExpr>(re->index)) {

        // if index constant, then add to set of constraints operating
        // on that array (actually, don't add constraint, just set index)
        llvm::BitVector &dis = elements[array];
        if (dis.empty())
          dis.resize(array->size);
        dis.set((unsigned)CE->getZExtValue(32));
      } else {
        // Remove information that array has been marked partially
        auto it2 = elements.find(array);
        if (it2 != elements.end())
          elements.erase(it2);

        // Add array as fully symbolic
        wholeObjects.insert(array);
      }
    }
  }
  IndependentElementSet(const IndependentElementSet &ies) = delete;

  IndependentElementSet &operator=(const IndependentElementSet &ies) = delete;

  IndependentElementSet(IndependentElementSet &&set) = default;
  IndependentElementSet &operator=(IndependentElementSet &&set) = default;
  void print(llvm::raw_ostream &os) const {
    os << "{";
    bool first = true;
    for (const auto &array : wholeObjects) {
      if (first) {
        first = false;
      } else {
        os << ", ";
      }
      os << "MO" << array->name;
    }
    for (const auto &element : elements) {
      const auto &array = element.first;
      const auto &dis = element.second;

      if (first) {
        first = false;
      } else {
        os << ", ";
      }

      os << "MO" << array->name << " : ";
      bool first = true;
      os << "{";
      for (auto it = dis.set_bits_begin(), ie = dis.set_bits_end(); it != ie;
           ++it) {
        if (first) {
          first = false;
        } else {
          os << ",";
        }
        os << *it;
      }
      os << "}";
    }
    os << "}";
  }

  // more efficient when this is the smaller set
  bool intersects(const IndependentElementSet &b) {
    // If there are any symbolic arrays in our query that b accesses
    for (const auto &array: wholeObjects) {
      if (b.wholeObjects.count(array) || 
          b.elements.find(array) != b.elements.end())
        return true;
    }
    for (auto &element : elements) {
      const Array *array = element.first;
      // if the array we access is symbolic in b
      if (b.wholeObjects.count(array))
        return true;

      auto it2 = b.elements.find(array);
      if (it2 == b.elements.end())
        continue;

      if (element.second.anyCommon(it2->second))
        return true;
    }
    return false;
  }

  // returns true iff set is changed by addition
  bool add(const IndependentElementSet &b) {
    exprs.insert(exprs.end(), b.exprs.begin(), b.exprs.end());

    bool modified = false;

    for (auto &array: b.wholeObjects) {
      elements_ty::iterator it2 = elements.find(array);
      if (it2!=elements.end()) {
        modified = true;
        elements.erase(it2);
        wholeObjects.insert(array);
      } else {
        if (!wholeObjects.count(array)) {
          modified = true;
          wholeObjects.insert(array);
        }
      }
    }
    for (auto &element : b.elements) {
      const auto &array = element.first;

      // Skip, if already a fully symbolic array
      if (wholeObjects.count(array))
        continue;

      // Try to insert element
      auto inserted = elements.insert(element);
      if (inserted.second)
        // Success
        modified = true;
      else {
        // Element already exists
        if (inserted.first->second.anyCommon(element.second)) {
          inserted.first->second |= element.second;
          modified = true;
        }
      }
    }
    return modified;
  }
};

inline llvm::raw_ostream &operator<<(llvm::raw_ostream &os,
                                     const IndependentElementSet &ies) {
  ies.print(os);
  return os;
}

// Breaks down a constraint into all of it's individual pieces, returning a
// list of IndependentElementSets or the independent factors.
std::vector<IndependentElementSet>
getAllIndependentConstraintsSets(const Query &query) {
  std::vector<IndependentElementSet> factors;

  ConstantExpr *CE = dyn_cast<ConstantExpr>(query.expr);
  if (CE) {
    assert(CE && CE->isFalse() && "the expr should always be false and "
                                  "therefore not included in factors");
  } else {
    ref<Expr> neg = Expr::createIsZero(query.expr);
    factors.emplace_back(IndependentElementSet(neg));
  }

  // iterate through all the previously separated constraints.  Until we
  // actually return, factors is treated as a queue of expressions to be
  // evaluated.  If the queue property isn't maintained, then the exprs
  // could be returned in an order different from how they came it, negatively
  // affecting later stages.
  for (const auto &constraint: query.constraints)
    factors.emplace_back(IndependentElementSet(constraint));

  bool doneLoop = false;
  do {
    doneLoop = true;
    std::vector<IndependentElementSet> done;
    while (!factors.empty()) {
      IndependentElementSet current = std::move(factors.back());
      factors.pop_back();
      // This list represents the set of factors that are separate from current.
      // Those that are not inserted into this list (queue) intersect with
      // current.
      std::vector<IndependentElementSet> keep;

      while (!factors.empty()) {
        IndependentElementSet compare = std::move(factors.back());
        factors.pop_back();
        if (current.intersects(compare)) {
          if (current.add(compare)) {
            // Means that we have added (z=y)added to (x=y)
            // Now need to see if there are any (z=?)'s
            doneLoop = false;
          }
        } else {
          keep.emplace_back(std::move(compare));
        }
      }
      done.emplace_back(std::move(current));
      factors.swap(keep);
    }
    factors = std::move(done);
  } while (!doneLoop);
  return factors;
}

static 
IndependentElementSet getIndependentConstraints(const Query& query,
                                                std::vector< ref<Expr> > &result) {
  IndependentElementSet eltsClosure(query.expr);

  std::vector< std::pair<ref<Expr>, IndependentElementSet> > worklist;
  for (auto &constraint : query.constraints)
    worklist.emplace_back(std::make_pair(constraint, IndependentElementSet(constraint)));

  // XXX This should be more efficient (in terms of low level copy stuff).
  bool done = false;
  do {
    done = true;
    std::vector< std::pair<ref<Expr>, IndependentElementSet> > newWorklist;
    for (auto &item: worklist) {
      if (item.second.intersects(eltsClosure)) {
        if (eltsClosure.add(item.second))
          done = false;
        result.push_back(item.first);
        // Means that we have added (z=y)added to (x=y)
        // Now need to see if there are any (z=?)'s
      } else {
        newWorklist.emplace_back(std::move(item));
      }
    }
    worklist.swap(newWorklist);
  } while (!done);

  KLEE_DEBUG(
    std::set< ref<Expr> > reqset(result.begin(), result.end());
    errs() << "--\n";
    errs() << "Q: " << query.expr << "\n";
    errs() << "\telts: " << IndependentElementSet(query.expr) << "\n";
    int i = 0;
    for (ConstraintManager::const_iterator it = query.constraints.begin(),
        ie = query.constraints.end(); it != ie; ++it) {
      errs() << "C" << i++ << ": " << *it;
      errs() << " " << (reqset.count(*it) ? "(required)" : "(independent)") << "\n";
      errs() << "\telts: " << IndependentElementSet(*it) << "\n";
    }
    errs() << "elts closure: " << eltsClosure << "\n";
 );


  return eltsClosure;
}


// Extracts which arrays are referenced from a particular independent set.  Examines both
// the actual known array accesses arr[1] plus the undetermined accesses arr[x].
static
void calculateArrayReferences(const IndependentElementSet & ie,
                              std::vector<const Array *> &returnVector){
  std::set<const Array*> thisSeen;
  for(auto &element: ie.elements){
    thisSeen.insert(element.first);
  }
  for(auto &element: ie.wholeObjects)
    thisSeen.insert(element);

  returnVector.insert(returnVector.end(), thisSeen.begin(), thisSeen.end());
}

class IndependentSolver : public SolverImpl {
private:
  Solver *solver;

public:
  IndependentSolver(Solver *_solver) 
    : solver(_solver) {}
  ~IndependentSolver() { delete solver; }

  bool computeTruth(const Query&, bool &isValid);
  bool computeValidity(const Query&, Solver::Validity &result);
  bool computeValue(const Query&, ref<Expr> &result);
  bool computeInitialValues(const Query& query,
                            const std::vector<const Array*> &objects,
                            std::vector< std::vector<unsigned char> > &values,
                            bool &hasSolution);
  SolverRunStatus getOperationStatusCode();
  char *getConstraintLog(const Query&);
  void setCoreSolverTimeout(time::Span timeout);
};
  
bool IndependentSolver::computeValidity(const Query& query,
                                        Solver::Validity &result) {
  std::vector< ref<Expr> > required;
  IndependentElementSet eltsClosure =
    getIndependentConstraints(query, required);
  ConstraintManager tmp(required);
  return solver->impl->computeValidity(Query(tmp, query.expr), 
                                       result);
}

bool IndependentSolver::computeTruth(const Query& query, bool &isValid) {
  std::vector< ref<Expr> > required;
  IndependentElementSet eltsClosure = 
    getIndependentConstraints(query, required);
  ConstraintManager tmp(required);
  return solver->impl->computeTruth(Query(tmp, query.expr), 
                                    isValid);
}

bool IndependentSolver::computeValue(const Query& query, ref<Expr> &result) {
  std::vector< ref<Expr> > required;
  IndependentElementSet eltsClosure = 
    getIndependentConstraints(query, required);
  ConstraintManager tmp(required);
  return solver->impl->computeValue(Query(tmp, query.expr), result);
}

// Helper function used only for assertions to make sure point created
// during computeInitialValues is in fact correct. The ``retMap`` is used
// in the case ``objects`` doesn't contain all the assignments needed.
bool assertCreatedPointEvaluatesToTrue(const Query &query,
                                       const std::vector<const Array*> &objects,
                                       std::vector< std::vector<unsigned char> > &values,
                                       std::map<const Array*, std::vector<unsigned char> > &retMap){
  // _allowFreeValues is set to true so that if there are missing bytes in the assigment
  // we will end up with a non ConstantExpr after evaluating the assignment and fail
  Assignment assign = Assignment(objects, values, /*_allowFreeValues=*/true);

  // Add any additional bindings.
  // The semantics of std::map should be to not insert a (key, value)
  // pair if it already exists so we should continue to use the assignment
  // from ``objects`` and ``values``.
  if (retMap.size() > 0)
    assign.bindings.insert(retMap.begin(), retMap.end());

  for(ConstraintManager::constraint_iterator it = query.constraints.begin();
      it != query.constraints.end(); ++it){
    ref<Expr> ret = assign.evaluate(*it);

    assert(isa<ConstantExpr>(ret) && "assignment evaluation did not result in constant");
    ref<ConstantExpr> evaluatedConstraint = dyn_cast<ConstantExpr>(ret);
    if(evaluatedConstraint->isFalse()){
      return false;
    }
  }
  ref<Expr> neg = Expr::createIsZero(query.expr);
  ref<Expr> q = assign.evaluate(neg);
  assert(isa<ConstantExpr>(q) && "assignment evaluation did not result in constant");
  return cast<ConstantExpr>(q)->isTrue();
}

bool IndependentSolver::computeInitialValues(const Query& query,
                                             const std::vector<const Array*> &objects,
                                             std::vector< std::vector<unsigned char> > &values,
                                             bool &hasSolution){
  // We assume the query has a solution except proven differently
  // This is important in case we don't have any constraints but
  // we need initial values for requested array objects.
  hasSolution = true;
  auto factors = getAllIndependentConstraintsSets(query);

  //Used to rearrange all of the answers into the correct order
  std::map<const Array*, std::vector<unsigned char> > retMap;
  for (auto &factor: factors) {
    std::vector<const Array*> arraysInFactor;
    calculateArrayReferences(factor, arraysInFactor);
    // Going to use this as the "fresh" expression for the Query() invocation below
    assert(factor.exprs.size() >= 1 && "No null/empty factors");
    if (arraysInFactor.empty())
      continue;

    ConstraintManager tmp(factor.exprs);
    std::vector<std::vector<unsigned char> > tempValues;
    if (!solver->impl->computeInitialValues(Query(tmp, ConstantExpr::alloc(0, Expr::Bool)),
                                            arraysInFactor, tempValues, hasSolution)){
      values.clear();
      return false;
    } else if (!hasSolution){
      values.clear();
      return true;
    } else {
      assert(tempValues.size() == arraysInFactor.size() &&
             "Should be equal number arrays and answers");
      for (unsigned i = 0; i < tempValues.size(); i++){
        if (retMap.count(arraysInFactor[i])){
          // We already have an array with some partially correct answers,
          // so we need to place the answers to the new query into the right
          // spot while avoiding the undetermined values also in the array
          std::vector<unsigned char> * tempPtr = &retMap[arraysInFactor[i]];
          assert(tempPtr->size() == tempValues[i].size() &&
                 "we're talking about the same array here");
          auto &ds = factor.elements[arraysInFactor[i]];
          for (auto it2 = ds.set_bits_begin(), it2E = ds.set_bits_end();
               it2 != it2E; it2++) {
            unsigned index = * it2;
            (* tempPtr)[index] = tempValues[i][index];
          }
        } else {
          // Dump all the new values into the array
          retMap[arraysInFactor[i]] = tempValues[i];
        }
      }
    }
  }
  for (std::vector<const Array *>::const_iterator it = objects.begin();
       it != objects.end(); it++){
    const Array * arr = * it;
    if (!retMap.count(arr)){
      // this means we have an array that is somehow related to the
      // constraint, but whose values aren't actually required to
      // satisfy the query.
      std::vector<unsigned char> ret(arr->size);
      values.push_back(ret);
    } else {
      values.push_back(retMap[arr]);
    }
  }
  assert(assertCreatedPointEvaluatesToTrue(query, objects, values, retMap) && "should satisfy the equation");
  return true;
}

SolverImpl::SolverRunStatus IndependentSolver::getOperationStatusCode() {
  return solver->impl->getOperationStatusCode();      
}

char *IndependentSolver::getConstraintLog(const Query& query) {
  return solver->impl->getConstraintLog(query);
}

void IndependentSolver::setCoreSolverTimeout(time::Span timeout) {
  solver->impl->setCoreSolverTimeout(timeout);
}

Solver *klee::createIndependentSolver(Solver *s) {
  return new Solver(new IndependentSolver(s));
}
