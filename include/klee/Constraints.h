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

// FIXME: Currently we use ConstraintManager for two things: to pass
// sets of constraints around, and to optimize constraints. We should
// move the first usage into a separate data structure
// (ConstraintSet?) which ConstraintManager could embed if it likes.
namespace klee {

class ExprVisitor;
  
class ConstraintManager {
public:
  typedef std::vector< ref<Expr> > constraints_ty;
  typedef constraints_ty::iterator iterator;
  typedef constraints_ty::const_iterator const_iterator;

  ConstraintManager() {}

  // create from constraints with no optimization
  explicit ConstraintManager(const std::vector<ref<Expr> >& _constraints);

  typedef std::vector< ref<Expr> >::const_iterator constraint_iterator;

  ref<Expr> simplifyExpr(ref<Expr> e) const;

  void addConstraint(ref<Expr> e);
  
  bool empty() const;
  ref<Expr> back() const;
  constraint_iterator begin() const;
  constraint_iterator end() const;
  size_t size() const;

  bool operator ==(const ConstraintManager& other) const;
  
private:
  std::vector< ref<Expr> > constraints;

  // returns true iff the constraints were modified
  bool rewriteConstraints(ExprVisitor &visitor);

  void addConstraintInternal(ref<Expr> e);
};

}

#endif /* KLEE_CONSTRAINTS_H */
