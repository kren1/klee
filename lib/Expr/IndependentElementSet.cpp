/*
 * IndependentElementSet.cpp
 *
 *  Created on: 30 Jan 2019
 *      Author: martin
 */

#include "klee/IndependentElementSet.h"

#include "klee/util/ExprUtil.h"

using namespace klee;

llvm::raw_ostream &klee::operator<<(llvm::raw_ostream &os,
                                    const IndependentElementSet &ies) {
  ies.print(os);
  return os;
}

IndependentElementSet::IndependentElementSet(ref<Expr> e) {
  exprs.push_back(e);
  // Track all reads in the program.  Determines whether reads are
  // concrete or symbolic.  If they are symbolic, "collapses" array
  // by adding it to wholeObjects.  Otherwise, creates a mapping of
  // the form Map<array, set<index>> which tracks which parts of the
  // array are being accessed.
  std::vector<ref<ReadExpr>> reads;
  findReads(e, /* visitUpdates= */
            true, reads);
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

      dis.set((unsigned)(CE->getZExtValue(32)));
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

void IndependentElementSet::print(llvm::raw_ostream &os) const {
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

bool IndependentElementSet::intersects(const IndependentElementSet &b) {
  // If there are any symbolic arrays in our query that b accesses
  for (const auto &array : wholeObjects) {
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

bool IndependentElementSet::add(const IndependentElementSet &b) {
  exprs.insert(exprs.end(), b.exprs.begin(), b.exprs.end());
  bool modified = false;
  for (auto &array : b.wholeObjects) {
    elements_ty::iterator it2 = elements.find(array);
    if (it2 != elements.end()) {
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
