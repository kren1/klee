//===-- ExprUtil.h ----------------------------------------------*- C++ -*-===//
//
//                     The KLEE Symbolic Virtual Machine
//
// This file is distributed under the University of Illinois Open Source
// License. See LICENSE.TXT for details.
//
//===----------------------------------------------------------------------===//

#ifndef KLEE_EXPRUTIL_H
#define KLEE_EXPRUTIL_H

#include "klee/Expr/ExprVisitor.h"
#include "klee/Expr/ExprVisitorT.h"

#include <vector>

namespace klee {
  class Array;
  class Expr;
  class ReadExpr;
  template<typename T> class ref;

  /// Find all ReadExprs used in the expression DAG. If visitUpdates
  /// is true then this will including those reachable by traversing
  /// update lists. Note that this may be slow and return a large
  /// number of results.
  void findReads(ref<Expr> e, 
                 bool visitUpdates,
                 std::vector< ref<ReadExpr> > &result);
  
  void findReads(ref<Expr> e, 
                 bool visitUpdates,
                 std::vector< ref<ReadExpr> > &result,
                 ExprHashSet &visited,
                 std::set<const UpdateNode *> &updates
                 );
  /// Return a list of all unique symbolic objects referenced by the given
  /// expression.
  void findSymbolicObjects(ref<Expr> e,
                           std::vector<const Array*> &results);

  /// Return a list of all unique symbolic objects referenced by the
  /// given expression range.
  template<typename InputIterator>
  void findSymbolicObjects(InputIterator begin, 
                           InputIterator end,
                           std::vector<const Array*> &results);

  class ConstantArrayFinderD : public ExprVisitor {
  protected:
    ExprVisitor::Action visitRead(const ReadExpr &re);

  public:
    std::set<const Array *> results;
  };

  class ConstantArrayFinderT : public ExprVisitorT<ConstantArrayFinderT> {
    friend ExprVisitorT<ConstantArrayFinderT>;
    typedef typename ExprVisitorT<ConstantArrayFinderT>::ActionT ActionT;

  protected:
    ActionT visitRead(const ReadExpr &re);

  public:
    std::set<const Array *> results;
  };
  using ConstantArrayFinder = ConstantArrayFinderT;
}

#endif /* KLEE_EXPRUTIL_H */
