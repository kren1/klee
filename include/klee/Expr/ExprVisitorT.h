//===-- ExprVisitor.h -------------------------------------------*- C++ -*-===//
//
//                     The KLEE Symbolic Virtual Machine
//
// This file is distributed under the University of Illinois Open Source
// License. See LICENSE.TXT for details.
//
//===----------------------------------------------------------------------===//

#ifndef KLEE_EXPRVISITORT_H
#define KLEE_EXPRVISITORT_H

#include "ExprHashMap.h"
#include "llvm/Support/CommandLine.h"

namespace klee {
extern llvm::cl::opt<bool> UseVisitorHash;
template <class T> class ExprVisitorT {
protected:
  // typed variant, but non-virtual for efficiency
  class ActionT {
  public:
    enum Kind { SkipChildren, DoChildren, ChangeTo };

  private:
    ActionT(Kind _kind) : kind(_kind), argument(nullptr) {}
    ActionT(Kind _kind, const ref<Expr> &_argument)
        : kind(_kind), argument(_argument) {}

  public:
    Kind kind;
    ref<Expr> argument;

    static ActionT changeTo(const ref<Expr> &expr) {
      return ActionT(ChangeTo, expr);
    }
    static ActionT doChildren() { return ActionT(DoChildren); }
    static ActionT skipChildren() { return ActionT(SkipChildren); }
  };

private:
  T &derived() { return *static_cast<T *>(this); }
  ref<Expr> visitActual(const ref<Expr> &e) {
    if (isa<ConstantExpr>(e)) {
      return e;
    } else {
      Expr &ep = *e.get();

      ActionT res = derived().visitExpr(ep);
      switch (res.kind) {
      case ActionT::DoChildren:
        // continue with normal action
        break;
      case ActionT::SkipChildren:
        return e;
      case ActionT::ChangeTo:
        return res.argument;
      }

      switch (ep.getKind()) {
      case Expr::NotOptimized:
        res = derived().visitNotOptimized(static_cast<NotOptimizedExpr &>(ep));
        break;
      case Expr::Read:
        res = derived().visitRead(static_cast<ReadExpr &>(ep));
        break;
      case Expr::Select:
        res = derived().visitSelect(static_cast<SelectExpr &>(ep));
        break;
      case Expr::Concat:
        res = derived().visitConcat(static_cast<ConcatExpr &>(ep));
        break;
      case Expr::Extract:
        res = derived().visitExtract(static_cast<ExtractExpr &>(ep));
        break;
      case Expr::ZExt:
        res = derived().visitZExt(static_cast<ZExtExpr &>(ep));
        break;
      case Expr::SExt:
        res = derived().visitSExt(static_cast<SExtExpr &>(ep));
        break;
      case Expr::Add:
        res = derived().visitAdd(static_cast<AddExpr &>(ep));
        break;
      case Expr::Sub:
        res = derived().visitSub(static_cast<SubExpr &>(ep));
        break;
      case Expr::Mul:
        res = derived().visitMul(static_cast<MulExpr &>(ep));
        break;
      case Expr::UDiv:
        res = derived().visitUDiv(static_cast<UDivExpr &>(ep));
        break;
      case Expr::SDiv:
        res = derived().visitSDiv(static_cast<SDivExpr &>(ep));
        break;
      case Expr::URem:
        res = derived().visitURem(static_cast<URemExpr &>(ep));
        break;
      case Expr::SRem:
        res = derived().visitSRem(static_cast<SRemExpr &>(ep));
        break;
      case Expr::Not:
        res = derived().visitNot(static_cast<NotExpr &>(ep));
        break;
      case Expr::And:
        res = derived().visitAnd(static_cast<AndExpr &>(ep));
        break;
      case Expr::Or:
        res = derived().visitOr(static_cast<OrExpr &>(ep));
        break;
      case Expr::Xor:
        res = derived().visitXor(static_cast<XorExpr &>(ep));
        break;
      case Expr::Shl:
        res = derived().visitShl(static_cast<ShlExpr &>(ep));
        break;
      case Expr::LShr:
        res = derived().visitLShr(static_cast<LShrExpr &>(ep));
        break;
      case Expr::AShr:
        res = derived().visitAShr(static_cast<AShrExpr &>(ep));
        break;
      case Expr::Eq:
        res = derived().visitEq(static_cast<EqExpr &>(ep));
        break;
      case Expr::Ne:
        res = derived().visitNe(static_cast<NeExpr &>(ep));
        break;
      case Expr::Ult:
        res = derived().visitUlt(static_cast<UltExpr &>(ep));
        break;
      case Expr::Ule:
        res = derived().visitUle(static_cast<UleExpr &>(ep));
        break;
      case Expr::Ugt:
        res = derived().visitUgt(static_cast<UgtExpr &>(ep));
        break;
      case Expr::Uge:
        res = derived().visitUge(static_cast<UgeExpr &>(ep));
        break;
      case Expr::Slt:
        res = derived().visitSlt(static_cast<SltExpr &>(ep));
        break;
      case Expr::Sle:
        res = derived().visitSle(static_cast<SleExpr &>(ep));
        break;
      case Expr::Sgt:
        res = derived().visitSgt(static_cast<SgtExpr &>(ep));
        break;
      case Expr::Sge:
        res = derived().visitSge(static_cast<SgeExpr &>(ep));
        break;
      case Expr::Constant:
      default:
        assert(0 && "invalid expression kind");
      }

      switch (res.kind) {
      default:
        assert(0 && "invalid kind");
      case ActionT::DoChildren: {
        bool rebuild = false;
        ref<Expr> e(&ep), kids[8];
        unsigned count = ep.getNumKids();
        for (unsigned i = 0; i < count; i++) {
          ref<Expr> kid = ep.getKid(i);
          kids[i] = visit(kid);
          if (kids[i] != kid)
            rebuild = true;
        }
        if (rebuild) {
          e = ep.rebuild(kids);
          if (recursive)
            e = visit(e);
        }
        if (!isa<ConstantExpr>(e)) {
          res = derived().visitExprPost(*e.get());
          if (res.kind == ActionT::ChangeTo)
            e = res.argument;
        }
        return e;
      }
      case ActionT::SkipChildren:
        return e;
      case ActionT::ChangeTo:
        return res.argument;
      }
    }
  }

protected:
  ActionT visitExpr(const Expr &) { return ActionT::doChildren(); }

  ActionT visitExprPost(const Expr &) { return ActionT::skipChildren(); }

  ActionT visitNotOptimized(const NotOptimizedExpr &) {
    return ActionT::doChildren();
  }

  ActionT visitRead(const ReadExpr &) { return ActionT::doChildren(); }

  ActionT visitSelect(const SelectExpr &) { return ActionT::doChildren(); }

  ActionT visitConcat(const ConcatExpr &) { return ActionT::doChildren(); }

  ActionT visitExtract(const ExtractExpr &) { return ActionT::doChildren(); }

  ActionT visitZExt(const ZExtExpr &) { return ActionT::doChildren(); }

  ActionT visitSExt(const SExtExpr &) { return ActionT::doChildren(); }

  ActionT visitAdd(const AddExpr &) { return ActionT::doChildren(); }

  ActionT visitSub(const SubExpr &) { return ActionT::doChildren(); }

  ActionT visitMul(const MulExpr &) { return ActionT::doChildren(); }

  ActionT visitUDiv(const UDivExpr &) { return ActionT::doChildren(); }

  ActionT visitSDiv(const SDivExpr &) { return ActionT::doChildren(); }

  ActionT visitURem(const URemExpr &) { return ActionT::doChildren(); }

  ActionT visitSRem(const SRemExpr &) { return ActionT::doChildren(); }

  ActionT visitNot(const NotExpr &) { return ActionT::doChildren(); }

  ActionT visitAnd(const AndExpr &) { return ActionT::doChildren(); }

  ActionT visitOr(const OrExpr &) { return ActionT::doChildren(); }

  ActionT visitXor(const XorExpr &) { return ActionT::doChildren(); }

  ActionT visitShl(const ShlExpr &) { return ActionT::doChildren(); }

  ActionT visitLShr(const LShrExpr &) { return ActionT::doChildren(); }

  ActionT visitAShr(const AShrExpr &) { return ActionT::doChildren(); }

  ActionT visitEq(const EqExpr &) { return ActionT::doChildren(); }

  ActionT visitNe(const NeExpr &) { return ActionT::doChildren(); }

  ActionT visitUlt(const UltExpr &) { return ActionT::doChildren(); }

  ActionT visitUle(const UleExpr &) { return ActionT::doChildren(); }

  ActionT visitUgt(const UgtExpr &) { return ActionT::doChildren(); }

  ActionT visitUge(const UgeExpr &) { return ActionT::doChildren(); }

  ActionT visitSlt(const SltExpr &) { return ActionT::doChildren(); }

  ActionT visitSle(const SleExpr &) { return ActionT::doChildren(); }

  ActionT visitSgt(const SgtExpr &) { return ActionT::doChildren(); }

  ActionT visitSge(const SgeExpr &) { return ActionT::doChildren(); }
  explicit ExprVisitorT(bool _recursive = false) : recursive(_recursive) {
    useVisitorHash = UseVisitorHash;
  }
  ExprVisitorT(bool visitorHash, bool _recursive = false)
      : recursive(_recursive), useVisitorHash(visitorHash) {}
  virtual ~ExprVisitorT() {}

private:
  typedef ExprHashMap<ref<Expr>> visited_ty;
  visited_ty visited;
  bool recursive;
  bool useVisitorHash = false;

public:
  // apply the visitor to the expression and return a possibly
  // modified new expression.
  ref<Expr> visit(const ref<Expr> &e) {
    if (!useVisitorHash || isa<ConstantExpr>(e)) {
      return visitActual(e);
    } else {
      visited_ty::iterator it = visited.find(e);

      if (it != visited.end()) {
        return it->second;
      } else {
        ref<Expr> res = visitActual(e);
        visited.insert(std::make_pair(e, res));
        return res;
      }
    }
  }
};
}

#endif /* KLEE_EXPRVISITORT_H */
