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

#include <vector>

using namespace klee;

//Root node belongs to all sets, so all 3 bits are set.
PTree::PTree(const data_type &root) : root(PTreeNodePtr(new Node(nullptr, root),7 )) {}

std::pair<PTreeNode*, PTreeNode*>
PTree::split(Node *n, 
             const data_type &leftData, 
             const data_type &rightData) {
  assert(n && !n->left.getPointer() && !n->right.getPointer());
  n->left = PTreeNodePtr(new Node(n, leftData));
  n->right = PTreeNodePtr(new Node(n, rightData));
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
    delete n;
    n = p;
  } while (n && !n->left.getPointer() && !n->right.getPointer());
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

