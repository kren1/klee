/*
 * IndependentElementSet.h
 *
 *  Created on: 17 Jan 2019
 *      Author: martin
 */

#ifndef INCLUDE_KLEE_INDEPENDENTELEMENTSET_H_
#define INCLUDE_KLEE_INDEPENDENTELEMENTSET_H_

#include "klee/Expr.h"
#include "llvm/ADT/BitVector.h"
#include <map>

namespace klee {
class IndependentElementSet {
public:
  typedef std::map<const Array *, llvm::BitVector> elements_ty;
  elements_ty
      elements; // Represents individual elements of array accesses (arr[1])
  std::set<const Array *>
      wholeObjects; // Represents symbolically accessed arrays (arr[x])
  std::vector<ref<Expr>> exprs; // All expressions that are associated with this
                                // factor Although order doesn't matter, we use
                                // a vector to match the ConstraintManager
                                // constructor that will eventually be invoked.

  IndependentElementSet() = delete;
  IndependentElementSet(ref<Expr> e);
  IndependentElementSet(const IndependentElementSet &ies) = delete;

  IndependentElementSet &operator=(const IndependentElementSet &ies) = delete;

  IndependentElementSet(IndependentElementSet &&set) = default;
  IndependentElementSet &operator=(IndependentElementSet &&set) = default;

  void print(llvm::raw_ostream &os) const;

  // more efficient when this is the smaller set
  bool intersects(const IndependentElementSet &b);

  // returns true iff set is changed by addition
  bool add(const IndependentElementSet &b);
};

llvm::raw_ostream &operator<<(llvm::raw_ostream &os,
                              const IndependentElementSet &ies);

} // namespace klee

#endif /* INCLUDE_KLEE_INDEPENDENTELEMENTSET_H_ */
