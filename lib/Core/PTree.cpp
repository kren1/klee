//===-- PTree.cpp ---------------------------------------------------------===//
//
//                     The KLEE Symbolic Virtual Machine
//
// This file is distributed under the University of Illinois Open Source
// License. See LICENSE.TXT for details.
//
//===----------------------------------------------------------------------===//

#include "PTree.h"

#include "klee/Expr/Expr.h"
#include "klee/Expr/ExprPPrinter.h"
#include "klee/Statistics.h"

#include <vector>

using namespace klee;

PTree::PTree(const data_type &_root) : root(PTreeNodePtr(new Node(nullptr, _root),0 )) {
    
    _root->ptreeNode = root.getPointer();
}
SQLIntStatistic totalStates("TotalStates", "TStates");

std::pair<PTreeNode*, PTreeNode*>
PTree::split(Node *n, 
             const data_type &leftData, 
             const data_type &rightData) {
  ++totalStates;
  auto node = n;
  assert(n && !n->left.getPointer() && !n->right.getPointer());
  assert(n == rightData->ptreeNode && "Split expect current to be right");
  uint8_t currentNodeTag = root.getInt();
  if (node->parent)
       currentNodeTag = node->parent->left.getPointer() == node
                                ? node->parent->left.getInt()
                                : node->parent->right.getInt();
  n->left = PTreeNodePtr(new(alloc.Allocate()) Node(n, leftData));
  n->right = PTreeNodePtr(new(alloc.Allocate()) Node(n, rightData), currentNodeTag);
  assert(n->left.getPointer() != n);
  assert(n->right.getPointer() != n);
  return std::make_pair(n->left.getPointer(), n->right.getPointer());
}

void PTree::remove(Node *n) {
  assert(!n->left.getPointer() && !n->right.getPointer());
  do {
    Node *p = n->parent;
    if (p) {
      if (n == p->left.getPointer()) {
        p->left = PTreeNodePtr(nullptr);
      } else {
        assert(n == p->right.getPointer());
        p->right = PTreeNodePtr(nullptr);
      }
    }
//    delete n;
    n->~Node();
    alloc.Deallocate(n);
    n = p;
  } while (n && !n->left.getPointer() && !n->right.getPointer());

  if (n) {
    // We're now at a node that has exactly one child; we've just deleted the
    // other one. Eliminate the node and connect its child to the parent
    // directly (if it's not the root).
    auto child = n->left.getPointer() ? n->left : n->right;
    Node *parent = n->parent;

    child.getPointer()->parent = parent;
    if (!parent) {
      // We're at the root.
      root = child;
  //    root.setInt(7);
    } else {
      if (n == parent->left.getPointer()) {
        parent->left = child;
      } else {
        assert(n == parent->right.getPointer());
        parent->right = child;
      }
    }

    n->~Node();
    alloc.Deallocate(n);
//    delete n;
  }
}

void PTree::dump(llvm::raw_ostream &os) {
  ExprPPrinter *pp = ExprPPrinter::create(os);
  pp->setNewline("\\l");
  os << "digraph G {\n";
  os << "\tsize=\"10,7.5\";\n";
  os << "\tratio=fill;\n";
  os << "\trotate=90;\n";
  os << "\tcenter = \"true\";\n";
  os << "\tnode [style=\"filled\",width=.1,height=.1,fontname=\"Terminus\"]\n";
  os << "\tedge [arrowsize=.3]\n";
  std::vector<PTree::Node*> stack;
  stack.push_back(root.getPointer());
  while (!stack.empty()) {
    PTree::Node *n = stack.back();
    stack.pop_back();
    os << "\tn" << n << " [shape=diamond";
    if (n->data)
      os << ",fillcolor=green";
    os << "];\n";
    if (n->left.getPointer()) {
      os << "\tn" << n << " -> n" << n->left.getPointer();
      os << " [label=" << (int)n->left.getInt() << "];\n";
      stack.push_back(n->left.getPointer());
    }
    if (n->right.getPointer()) {
      os << "\tn" << n << " -> n" << n->right.getPointer();
      os << " [label=" << (int)n->right.getInt() << "];\n";
      stack.push_back(n->right.getPointer());
    }
  }
  os << "}\n";
  delete pp;
}

PTreeNode::PTreeNode(PTreeNode * parent, ExecutionState * data)
  : parent{parent}, data{data} {
    left = PTreeNodePtr(nullptr);
    right = PTreeNodePtr(nullptr);
  }

