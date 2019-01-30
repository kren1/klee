//===-- Constraints.h -------------------------------------------*- C++ -*-===//
//
//                     The KLEE Symbolic Virtual Machine
//
// This file is distributed under the University of Illinois Open Source
// License. See LICENSE.TXT for details.
//
//===----------------------------------------------------------------------===//

#ifndef KLEE_CONSTRAINTS_H
#define KLEE_CONSTRAINTS_H

#include "klee/Expr.h"

namespace klee {

/// Resembles a set of constraints that can be passed around
///
class ConstraintSet {
  friend class ConstraintManager;

public:
  using constraints_ty = std::vector<ref<Expr>>;
  typedef constraints_ty::iterator iterator;
  typedef constraints_ty::const_iterator const_iterator;
  typedef const_iterator constraint_iterator;

  virtual constraint_iterator begin() const = 0;
  virtual constraint_iterator end() const = 0;
  //  size_t size() const;

  ConstraintSet() = default;

  virtual ~ConstraintSet() = default;
};

class SimpleConstraintSet : public ConstraintSet {
public:
  virtual constraint_iterator begin() const { return constraints.begin(); }
  virtual constraint_iterator end() const { return constraints.end(); }

  bool operator==(const SimpleConstraintSet &b) const {
    return constraints == b.constraints;
  }

  void push_back(const ref<Expr> &e) { constraints.push_back(e); }

  SimpleConstraintSet(constraints_ty cs) : constraints(cs) {}
  SimpleConstraintSet(const ConstraintSet &cs) {
    for (auto &c : cs)
      constraints.push_back(c);
  }

  SimpleConstraintSet() = default;

protected:
  constraints_ty constraints;
};

class OrderedConstraintSet : public ConstraintSet {
public:
  virtual constraint_iterator begin() const { return constraints.begin(); }
  virtual constraint_iterator end() const { return constraints.end(); }

  bool operator==(const OrderedConstraintSet &b) const {
    return constraints == b.constraints;
  }

  OrderedConstraintSet(constraints_ty cs) : constraints(cs) {
    std::sort(constraints.begin(), constraints.end(), cmp);
  }

  OrderedConstraintSet() = default;

protected:
  void push_back(const ref<Expr> &item) {
    // Search for the correct position
    auto position =
        std::upper_bound(constraints.begin(), constraints.end(), item, cmp);
    constraints.insert(position, item);
  }

  /// Comparator for two expressions.
  ///
  /// @param a expression
  /// @param b expression
  /// @return true if expression a is semantically ordered before b, otherwise
  /// false
  static bool cmp(const ref<Expr> &a, const ref<Expr> &b) { return *a < *b; }

  constraints_ty constraints;
};

///// Manages constraints in sets of independent constraints
/////
// class IndepConstraintSet {
// public:
//  std::vector<IndependentElementSet> idep_constraints;
//};

class ExprVisitor;

/// Manages constraints, e.g. optimisation
class ConstraintManager {
public:
  // create from constraints with no optimization
  explicit ConstraintManager(SimpleConstraintSet &constraints);

  static ref<Expr> simplifyExpr(const ConstraintSet &constraints,
                                const ref<Expr> &e);

  void addConstraint(ref<Expr> e);

private:
  // returns true iff the constraints were modified
  bool rewriteConstraints(ExprVisitor &visitor);

  void addConstraintInternal(ref<Expr> e);

  SimpleConstraintSet &constraints;
};

}

#endif /* KLEE_CONSTRAINTS_H */
