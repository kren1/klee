/*
 * IndependentElementSet.h
 *
 *  Created on: 17 Jan 2019
 *      Author: martin
 */

#ifndef LIB_EXPR_INDEPENDENTELEMENTSET_H_
#define LIB_EXPR_INDEPENDENTELEMENTSET_H_

#include "klee/Expr.h"
#include "llvm/ADT/BitVector.h"
#include <map>

namespace klee {
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

}

#endif /* LIB_EXPR_INDEPENDENTELEMENTSET_H_ */
