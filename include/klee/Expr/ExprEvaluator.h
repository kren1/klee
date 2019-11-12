//===-- ExprEvaluator.h -----------------------------------------*- C++ -*-===//
//
//                     The KLEE Symbolic Virtual Machine
//
// This file is distributed under the University of Illinois Open Source
// License. See LICENSE.TXT for details.
//
//===----------------------------------------------------------------------===//

#ifndef KLEE_EXPREVALUATOR_H
#define KLEE_EXPREVALUATOR_H

#include "klee/Expr/Expr.h"
#include "klee/Expr/ExprVisitor.h"
#include "klee/Expr/ExprVisitorT.h"

namespace klee {
  class ExprEvaluator : public ExprVisitor {
  protected:
    Action evalRead(const UpdateList &ul, unsigned index);
    Action visitRead(const ReadExpr &re);
    Action visitExpr(const Expr &e);
      
    Action protectedDivOperation(const BinaryExpr &e);
    Action visitUDiv(const UDivExpr &e);
    Action visitSDiv(const SDivExpr &e);
    Action visitURem(const URemExpr &e);
    Action visitSRem(const SRemExpr &e);
    Action visitExprPost(const Expr& e);
      
  public:
    ExprEvaluator() {}

    /// getInitialValue - Return the initial value for a symbolic byte.
    ///
    /// This will only be called for constant arrays if the index is
    /// out-of-bounds. If the value is unknown then the user should return a
    /// ReadExpr at the initial version of this array.
    virtual ref<Expr> getInitialValue(const Array& os, unsigned index) = 0;
  };

  template <class T>
  class ExprEvaluatorT : public ExprVisitorT<ExprEvaluatorT<T>> {
    friend ExprVisitorT<ExprEvaluatorT<T>>;
    typedef typename ExprVisitorT<ExprEvaluatorT<T>>::ActionT ActionT;
    T &derived() { return *static_cast<T *>(this); }

  protected:
    ActionT evalRead(const UpdateList &ul, unsigned index) {
      for (const UpdateNode *un = ul.head; un; un = un->next) {
        ref<Expr> ui = derived().visit(un->index);

        if (ConstantExpr *CE = dyn_cast<ConstantExpr>(ui)) {
          if (CE->getZExtValue() == index)
            return ActionT::changeTo(derived().visit(un->value));
        } else {
          // update index is unknown, so may or may not be index, we
          // cannot guarantee value. we can rewrite to read at this
          // version though (mostly for debugging).

          return ActionT::changeTo(ReadExpr::create(
              UpdateList(ul.root, un),
              ConstantExpr::alloc(index, ul.root->getDomain())));
        }
      }

      if (ul.root->isConstantArray() && index < ul.root->size)
        return ActionT::changeTo(ul.root->constantValues[index]);

      return ActionT::changeTo(derived().getInitialValue(*ul.root, index));
    }

    ActionT visitExpr(const Expr &e) {
      // Evaluate all constant expressions here, in case they weren't folded in
      // construction. Don't do this for reads though, because we want them to
      // go to
      // the normal rewrite path.
      unsigned N = e.getNumKids();
      if (!N || isa<ReadExpr>(e))
        return ActionT::doChildren();

      for (unsigned i = 0; i != N; ++i)
        if (!isa<ConstantExpr>(e.getKid(i)))
          return ActionT::doChildren();

      ref<Expr> Kids[3];
      for (unsigned i = 0; i != N; ++i) {
        assert(i < 3);
        Kids[i] = e.getKid(i);
      }

      return ActionT::changeTo(e.rebuild(Kids));
    }

    ActionT visitRead(const ReadExpr &re) {
      ref<Expr> v = derived().visit(re.index);

      if (ConstantExpr *CE = dyn_cast<ConstantExpr>(v)) {
        return evalRead(re.updates, CE->getZExtValue());
      } else {
        return ActionT::doChildren();
      }
    }

    // we need to check for div by zero during partial evaluation,
    // if this occurs then simply ignore the 0 divisor and use the
    // original expression.
    ActionT protectedDivOperation(const BinaryExpr &e) {
      ref<Expr> kids[2] = {derived().visit(e.left), derived().visit(e.right)};

      if (ConstantExpr *CE = dyn_cast<ConstantExpr>(kids[1]))
        if (CE->isZero())
          kids[1] = e.right;

      if (kids[0] != e.left || kids[1] != e.right) {
        return ActionT::changeTo(e.rebuild(kids));
      } else {
        return ActionT::skipChildren();
      }
    }

    ActionT visitUDiv(const UDivExpr &e) { return protectedDivOperation(e); }
    ActionT visitSDiv(const SDivExpr &e) { return protectedDivOperation(e); }
    ActionT visitURem(const URemExpr &e) { return protectedDivOperation(e); }
    ActionT visitSRem(const SRemExpr &e) { return protectedDivOperation(e); }

    ActionT visitExprPost(const Expr &e) {
      // When evaluating an assignment we should fold NotOptimizedExpr
      // nodes so we can fully evaluate.
      if (e.getKind() == Expr::NotOptimized) {
        return ActionT::changeTo(static_cast<const NotOptimizedExpr &>(e).src);
      }
      return ActionT::skipChildren();
    }

  public:
    ExprEvaluatorT() {}
  };
}

#endif /* KLEE_EXPREVALUATOR_H */
