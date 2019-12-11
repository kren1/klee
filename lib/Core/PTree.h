//===-- PTree.h -------------------------------------------------*- C++ -*-===//
//
//                     The KLEE Symbolic Virtual Machine
//
// This file is distributed under the University of Illinois Open Source
// License. See LICENSE.TXT for details.
//
//===----------------------------------------------------------------------===//

#ifndef KLEE_PTREE_H
#define KLEE_PTREE_H

#include "klee/Expr/Expr.h"
#include "klee/ExecutionState.h"
#include "llvm/ADT/PointerIntPair.h"
#include "llvm/Support/RecyclingAllocator.h"
#include "llvm/Support/Allocator.h"

namespace klee {
  using PTreeNodePtr = llvm::PointerIntPair<PTreeNode*,3,uint8_t>;

  class PTreeNode {
    friend class PTree;
  public:
    PTreeNode *parent = nullptr;
    PTreeNodePtr left; 
    PTreeNodePtr right ;
    ExecutionState *data = nullptr;

  private:
    PTreeNode(PTreeNode * parent, ExecutionState * data);
    ~PTreeNode() = default;
  };

  using PTreeAlloc = llvm::RecyclingAllocator<llvm::BumpPtrAllocatorImpl<llvm::MallocAllocator,32768>, PTreeNode>;

  class PTree { 
    typedef ExecutionState* data_type;

  public:
    typedef class PTreeNode Node;
    PTreeNodePtr root;
    PTreeAlloc alloc;

    explicit PTree(const data_type &_root);
    ~PTree() = default;
    
    std::pair<Node*,Node*> split(Node *n,
                                 const data_type &leftData,
                                 const data_type &rightData);
    void remove(Node *n);

    void dump(llvm::raw_ostream &os);
  };


}

#endif /* KLEE_PTREE_H */
