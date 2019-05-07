//===-- Z3IntBuilder.cpp ------------------------------------------*- C++ -*-====//
//
//                     The KLEE Symbolic Virtual Machine
//
// This file is distributed under the University of Illinois Open Source
// License. See LICENSE.TXT for details.
//
//===----------------------------------------------------------------------===//
#include "klee/Config/config.h"
#ifdef ENABLE_Z3
#include "Z3IntBuilder.h"

#include "klee/Expr.h"
#include "klee/Internal/Support/ErrorHandling.h"
#include "klee/Solver.h"
#include "klee/SolverStats.h"
#include "klee/util/Bits.h"

#include "llvm/ADT/StringExtras.h"
#include "llvm/Support/CommandLine.h"

using namespace klee;

namespace {
llvm::cl::opt<bool> UseConstructHashZ3(
    "use-construct-hash-z3-int",
    llvm::cl::desc("Use hash-consing during Z3 query construction."),
    llvm::cl::init(true));

// FIXME: This should be std::atomic<bool>. Need C++11 for that.
bool Z3InterationLogOpen = false;
}

namespace klee {

// Declared here rather than `Z3IntBuilder.h` so they can be called in gdb.
template <> void Z3NodeHandle<Z3_sort>::dump();
template <> void Z3NodeHandle<Z3_ast>::dump();

void custom_z3_error_handler(Z3_context ctx, Z3_error_code ec);

Z3IntBuilder::Z3IntBuilder(bool autoClearConstructCache, const char* z3LogInteractionFileArg)
    : autoClearConstructCache(autoClearConstructCache), z3LogInteractionFile("") {
  if (z3LogInteractionFileArg)
    this->z3LogInteractionFile = std::string(z3LogInteractionFileArg);
  if (z3LogInteractionFile.length() > 0) {
    klee_message("Logging Z3 API interaction to \"%s\"",
                 z3LogInteractionFile.c_str());
    assert(!Z3InterationLogOpen && "interaction log should not already be open");
    Z3_open_log(z3LogInteractionFile.c_str());
    Z3InterationLogOpen = true;
  }
  // FIXME: Should probably let the client pass in a Z3_config instead
  Z3_config cfg = Z3_mk_config();
  // It is very important that we ask Z3 to let us manage memory so that
  // we are able to cache expressions and sorts.
  ctx = Z3_mk_context_rc(cfg);
  // Make sure we handle any errors reported by Z3.
  Z3_set_error_handler(ctx, custom_z3_error_handler);
  // When emitting Z3 expressions make them SMT-LIBv2 compliant
  Z3_set_ast_print_mode(ctx, Z3_PRINT_SMTLIB2_COMPLIANT);
  Z3_del_config(cfg);
}

Z3IntBuilder::~Z3IntBuilder() {
  // Clear caches so exprs/sorts gets freed before the destroying context
  // they aren associated with.
  clearConstructCache();
  _arr_hash.clear();
  constant_array_assertions.clear();
  Z3_del_context(ctx);
  if (z3LogInteractionFile.length() > 0) {
    Z3_close_log();
    Z3InterationLogOpen = false;
  }
}

Z3SortHandle Z3IntBuilder::getIntSort() {
  // FIXME: cache these
  return Z3SortHandle(Z3_mk_int_sort(ctx), ctx);
}

Z3SortHandle Z3IntBuilder::getArraySort(Z3SortHandle domainSort,
                                     Z3SortHandle rangeSort) {
  // FIXME: cache these
  return Z3SortHandle(Z3_mk_array_sort(ctx, domainSort, rangeSort), ctx);
}

Z3ASTHandle Z3IntBuilder::buildArray(const char *name, unsigned indexWidth,
                                  unsigned valueWidth) {
  Z3SortHandle domainSort = getIntSort();
  Z3SortHandle rangeSort = getIntSort();
  Z3SortHandle t = getArraySort(domainSort, rangeSort);
  Z3_symbol s = Z3_mk_string_symbol(ctx, const_cast<char *>(name));
  return Z3ASTHandle(Z3_mk_const(ctx, s, t), ctx);
}

Z3ASTHandle Z3IntBuilder::getTrue() { return Z3ASTHandle(Z3_mk_true(ctx), ctx); }

Z3ASTHandle Z3IntBuilder::getFalse() { return Z3ASTHandle(Z3_mk_false(ctx), ctx); }

Z3ASTHandle Z3IntBuilder::uIntConst(uint64_t value) {
  Z3SortHandle t = getIntSort();
  return Z3ASTHandle(Z3_mk_unsigned_int64(ctx, value, t), ctx);
}

Z3ASTHandle Z3IntBuilder::sIntConst(int64_t value) {
  Z3SortHandle t = getIntSort();
  return Z3ASTHandle(Z3_mk_int64(ctx, value, t), ctx);
}

Z3ASTHandle Z3IntBuilder::eqExpr(Z3ASTHandle a, Z3ASTHandle b) {
  return Z3ASTHandle(Z3_mk_eq(ctx, a, b), ctx);
}

Z3ASTHandle Z3IntBuilder::geExpr(Z3ASTHandle a, Z3ASTHandle b) {
  return Z3ASTHandle(Z3_mk_ge(ctx, a, b), ctx);
}

Z3ASTHandle Z3IntBuilder::ltExpr(Z3ASTHandle a, Z3ASTHandle b) {
  return Z3ASTHandle(Z3_mk_lt(ctx, a, b), ctx);
}

Z3ASTHandle Z3IntBuilder::leExpr(Z3ASTHandle a, Z3ASTHandle b) {
  return Z3ASTHandle(Z3_mk_le(ctx, a, b), ctx);
}
Z3ASTHandle Z3IntBuilder::notExpr(Z3ASTHandle expr) {
  return Z3ASTHandle(Z3_mk_not(ctx, expr), ctx);
}

Z3ASTHandle Z3IntBuilder::andExpr(Z3ASTHandle lhs, Z3ASTHandle rhs) {
  ::Z3_ast args[2] = {lhs, rhs};
  return Z3ASTHandle(Z3_mk_and(ctx, 2, args), ctx);
}

Z3ASTHandle Z3IntBuilder::orExpr(Z3ASTHandle lhs, Z3ASTHandle rhs) {
  ::Z3_ast args[2] = {lhs, rhs};
  return Z3ASTHandle(Z3_mk_or(ctx, 2, args), ctx);
}


Z3ASTHandle Z3IntBuilder::iffExpr(Z3ASTHandle lhs, Z3ASTHandle rhs) {
  Z3SortHandle lhsSort = Z3SortHandle(Z3_get_sort(ctx, lhs), ctx);
  Z3SortHandle rhsSort = Z3SortHandle(Z3_get_sort(ctx, rhs), ctx);
  assert(Z3_get_sort_kind(ctx, lhsSort) == Z3_get_sort_kind(ctx, rhsSort) &&
         "lhs and rhs sorts must match");
  assert(Z3_get_sort_kind(ctx, lhsSort) == Z3_BOOL_SORT && "args must have BOOL sort");
  return Z3ASTHandle(Z3_mk_iff(ctx, lhs, rhs), ctx);
}

Z3ASTHandle Z3IntBuilder::writeExpr(Z3ASTHandle array, Z3ASTHandle index,
                                 Z3ASTHandle value) {
  return Z3ASTHandle(Z3_mk_store(ctx, array, index, value), ctx);
}

Z3ASTHandle Z3IntBuilder::readExpr(Z3ASTHandle array, Z3ASTHandle index) {
  return Z3ASTHandle(Z3_mk_select(ctx, array, index), ctx);
}

Z3ASTHandle Z3IntBuilder::iteExpr(Z3ASTHandle condition, Z3ASTHandle whenTrue,
                               Z3ASTHandle whenFalse) {
  return Z3ASTHandle(Z3_mk_ite(ctx, condition, whenTrue, whenFalse), ctx);
}

Z3ASTHandle Z3IntBuilder::getInitialArray(const Array *root) {

  assert(root);
  Z3ASTHandle array_expr;
  bool hashed = _arr_hash.lookupArrayExpr(root, array_expr);

  if (!hashed) {
    // Unique arrays by name, so we make sure the name is unique by
    // using the size of the array hash as a counter.
    std::string unique_id = llvm::itostr(_arr_hash._array_hash.size());
    unsigned const uid_length = unique_id.length();
    unsigned const space = (root->name.length() > 32 - uid_length)
                               ? (32 - uid_length)
                               : root->name.length();
    std::string unique_name = root->name.substr(0, space) + unique_id;

    array_expr = buildArray(unique_name.c_str(), root->getDomain(),
                            root->getRange());

    if (root->isConstantArray() && constant_array_assertions.count(root) == 0) {
      std::vector<Z3ASTHandle> array_assertions;
      assert(root->valueType != Expr::InvalidWidth && "Int solver can't deal with unknown type constant arrays");
      unsigned byteStride = root->valueType / 8;
      for (unsigned i = 0, e = (root->size / byteStride); i != e; ++i) {
        // construct(= (select i root) root->value[i]) to be asserted in
        // Z3Solver.cpp
        int width_out;
	      
	// start with the MSB of the the stride
    ref<ConstantExpr> init = root->constantValues[((i + 1) * byteStride) - 1];
    for (size_t j = 1; j < byteStride; j++){
      // Concat does essentially a  lhs <<=| rhs operation
      ref<ConstantExpr> rhs = root->constantValues[((i + 1) * byteStride) - j - 1];
      init = init->Concat(rhs);
    }


        Z3ASTHandle array_value = construct(init, &width_out);
        auto rExpr = readExpr(array_expr, uIntConst(i)); 

        array_assertions.push_back(eqExpr(rExpr, array_value));
      }
      constant_array_assertions[root] = std::move(array_assertions);

    } else if (!root->isConstantArray()) {
      assert(root->valueType != Expr::InvalidWidth && "Int solver can't deal with unknown type symbolic arrays");
      std::vector<Z3ASTHandle> array_assertions;

      unsigned byteStride = root->valueType / 8;
//      llvm::errs() << root->name << ":" << root->valueType << "\n";
      for (unsigned i = 0, e = (root->size / byteStride); i != e; ++i) {
        auto rExpr = readExpr(array_expr, uIntConst(i)); 
        switch(byteStride) {
            case 1: 
                  array_assertions.push_back(leExpr(rExpr,sIntConst(std::numeric_limits<std::int8_t>::max())));
                  array_assertions.push_back(geExpr(rExpr,sIntConst(std::numeric_limits<std::int8_t>::min())));
                  break;
            case 2:
                  array_assertions.push_back(leExpr(rExpr,sIntConst(std::numeric_limits<std::int16_t>::max())));
                  array_assertions.push_back(geExpr(rExpr,sIntConst(std::numeric_limits<std::int16_t>::min())));
                  break;
            case 4: 
                  array_assertions.push_back(leExpr(rExpr,sIntConst(std::numeric_limits<std::int32_t>::max())));
                  array_assertions.push_back(geExpr(rExpr,sIntConst(std::numeric_limits<std::int32_t>::min())));
                  break;
            case 8:
                  array_assertions.push_back(leExpr(rExpr,sIntConst(std::numeric_limits<std::int64_t>::max())));
                  array_assertions.push_back(geExpr(rExpr,sIntConst(std::numeric_limits<std::int64_t>::min())));
                  break;
            default: assert(0 && "Invalid byte stride");
        }
        
      }
      constant_array_assertions[root] = std::move(array_assertions);
    }

    _arr_hash.hashArrayExpr(root, array_expr);
  }

  return (array_expr);
}

Z3ASTHandle Z3IntBuilder::getInitialRead(const Array *root, unsigned index) {
  return readExpr(getInitialArray(root), uIntConst(index));
}

Z3ASTHandle Z3IntBuilder::getArrayForUpdate(const Array *root,
                                         const UpdateNode *un) {
  if (!un) {
    return (getInitialArray(root));
  } else {
    // FIXME: This really needs to be non-recursive.
    Z3ASTHandle un_expr;
    bool hashed = _arr_hash.lookupUpdateNodeExpr(un, un_expr);

    if (!hashed) {
      un_expr = writeExpr(getArrayForUpdate(root, un->next),
                          construct(un->index, 0), construct(un->value, 0));

      _arr_hash.hashUpdateNodeExpr(un, un_expr);
    }

    return (un_expr);
  }
}

/** if *width_out!=1 then result is a bitvector,
    otherwise it is a bool */
Z3ASTHandle Z3IntBuilder::construct(ref<Expr> e, int *width_out) {
  // TODO: We could potentially use Z3_simplify() here
  // to store simpler expressions.
  if (!UseConstructHashZ3 || isa<ConstantExpr>(e)) {
    return constructActual(e, width_out);
  } else {
    ExprHashMap<std::pair<Z3ASTHandle, unsigned> >::iterator it =
        constructed.find(e);
    if (it != constructed.end()) {
      if (width_out)
        *width_out = it->second.second;
      return it->second.first;
    } else {
      int width;
      if (!width_out)
        width_out = &width;
      Z3ASTHandle res = constructActual(e, width_out);
      constructed.insert(std::make_pair(e, std::make_pair(res, *width_out)));
      return res;
    }
  }
}
bool isReadExprAtOffset(ref<Expr> e, const ReadExpr *base, ref<Expr> offset) {
  const ReadExpr *re = dyn_cast<ReadExpr>(e.get());
    
  // right now, all Reads are byte reads but some
  // transformations might change this
  if (!re || (re->getWidth() != Expr::Int8))
    return false;
    
  // Check if the index follows the stride. 
  // FIXME: How aggressive should this be simplified. The
  // canonicalizing builder is probably the right choice, but this
  // is yet another area where we would really prefer it to be
  // global or else use static methods.
  return SubExpr::create(re->index, base->index) == offset;
}
 
const ReadExpr* Z3IntBuilder::hasOrderedReads(const Expr *concatExpr, int stride, int& width) {
    assert(concatExpr->getKind() == Expr::Concat);
    assert(stride == 1 || stride == -1);
    
    const ReadExpr *base = dyn_cast<ReadExpr>(concatExpr->getKid(0));
    
    // right now, all Reads are byte reads but some
    // transformations might change this
    if (!base || base->getWidth() != Expr::Int8)
      return NULL;
    
    // Get stride expr in proper index width.
    width = 8;
    Expr::Width idxWidth = base->index->getWidth();
    ref<Expr> strideExpr = ConstantExpr::alloc(stride, idxWidth);
    ref<Expr> offset = ConstantExpr::create(0, idxWidth);
    
    width += 8;
    auto e = concatExpr->getKid(1);
    
    // concat chains are unbalanced to the right
    while (e->getKind() == Expr::Concat) {
      offset = AddExpr::create(offset, strideExpr);
      width += 8;
      if (!isReadExprAtOffset(e->getKid(0), base, offset))
	return NULL;
      
      e = e->getKid(1);
    }
    
    offset = AddExpr::create(offset, strideExpr);
    if (!isReadExprAtOffset(e, base, offset))
      return NULL;
    
    if (stride == -1) {
      auto re = cast<ReadExpr>(e.get());
      return re;
    }
    else return base;
}



class SignHandler {
    bool old;
    static bool sign_int;
public:
    SignHandler(bool val): old(sign_int){
        sign_int  = val;
    }
    ~SignHandler() {
        sign_int = old;
    }
    static bool isSign() { return sign_int; }
};
bool SignHandler::sign_int = false;
/** if *width_out!=1 then result is a bitvector,
    otherwise it is a bool */
Z3ASTHandle Z3IntBuilder::constructActual(ref<Expr> e, int *width_out) {
  int width;
  if (!width_out)
    width_out = &width;

  ++stats::queryConstructs;

  switch (e->getKind()) {
  case Expr::Constant: {
    ConstantExpr *CE = cast<ConstantExpr>(e);
    *width_out = CE->getWidth();

    // Coerce to bool if necessary.
    if (*width_out == 1)
      return CE->isTrue() ? getTrue() : getFalse();

    assert(*width_out <= 64 && "Unangled bigger than 64 bit wide constants");
    // Fast path.
    //Can't use the sign here, because the equality (and maybe other operations)will still be unsigned.
    return sIntConst(CE->getAPValue().getSExtValue());
    if (false && SignHandler::isSign())
      return sIntConst(CE->getAPValue().getSExtValue());
    else
      return uIntConst(CE->getZExtValue());

  }

  // Special
  case Expr::NotOptimized: {
    NotOptimizedExpr *noe = cast<NotOptimizedExpr>(e);
    return construct(noe->src, width_out);
  }

  case Expr::Read: {
    ReadExpr *re = cast<ReadExpr>(e);
    assert(re && re->updates.root);
    *width_out = re->updates.root->getRange();
    assert(re->getWidth() == re->updates.root->valueType && "Read of different width the the array type, meaning it is not an int read");
    return readExpr(getArrayForUpdate(re->updates.root, re->updates.head),
                    construct(re->index, 0));
  }

  case Expr::Select: {
    SelectExpr *se = cast<SelectExpr>(e);
    Z3ASTHandle cond = construct(se->cond, 0);
    Z3ASTHandle tExpr = construct(se->trueExpr, width_out);
    Z3ASTHandle fExpr = construct(se->falseExpr, width_out);
    return iteExpr(cond, tExpr, fExpr);
  }

  case Expr::Concat: {
    int readlsb_width;
    auto re = hasOrderedReads(e.get(), -1, readlsb_width);
    assert(re && "Int Solver can't handle non ordered reads");
    *width_out = readlsb_width;
    auto stride = re->updates.root->valueType / 8;
//    re->dump();
//    llvm::errs() << "of width: " << *width_out << " array type: " << re->updates.root->valueType << "\n";
    assert(readlsb_width == re->updates.root->valueType && "ReadLSB of different width the the array type, meaning it is not an int read");
    assert(stride > 0 && "can;t concatnact for unknwon array type");
    ref<Expr> index = UDivExpr::create(re->index, ConstantExpr::create(stride, re->index->getWidth()));
    return readExpr(getArrayForUpdate(re->updates.root, re->updates.head),
                    construct(index, 0));
 
  }

  case Expr::Extract: {
    assert(0 && "Can't handle extract");
  }

  // Casting
  case Expr::ZExt: 
  case Expr::SExt: {
    SignHandler sh(e->getKind() == Expr::SExt);
    CastExpr *ce = cast<CastExpr>(e);
    int srcWidth;
    Z3ASTHandle src = construct(ce->src, &srcWidth);
    *width_out = ce->getWidth();
    if(ce->src->getWidth() == Expr::Bool)  {
      return iteExpr(src, sIntConst(1), sIntConst(0));
    }
    return src;
  }

  // Arithmetic
  case Expr::Add: {
    AddExpr *ae = cast<AddExpr>(e);
    Z3ASTHandle left = construct(ae->left, width_out);
    Z3ASTHandle right = construct(ae->right, width_out);
    assert(*width_out != 1 && "uncanonicalized add");
    Z3_ast ar[2] = {left, right};
    Z3ASTHandle result = Z3ASTHandle(Z3_mk_add(ctx, 2, ar), ctx);
    return result;
  }

  case Expr::Sub: {
    SubExpr *se = cast<SubExpr>(e);
    Z3ASTHandle left = construct(se->left, width_out);
    Z3ASTHandle right = construct(se->right, width_out);
    assert(*width_out != 1 && "uncanonicalized sub");
    Z3_ast ar[2] = {left, right};
    Z3ASTHandle result = Z3ASTHandle(Z3_mk_sub(ctx, 2, ar), ctx);
    return result;
  }

  case Expr::Mul: {
    MulExpr *me = cast<MulExpr>(e);
    Z3ASTHandle right = construct(me->right, width_out);
    assert(*width_out != 1 && "uncanonicalized mul");
    Z3ASTHandle left = construct(me->left, width_out);
    Z3_ast ar[2] = {left, right};
    Z3ASTHandle result = Z3ASTHandle(Z3_mk_mul(ctx, 2, ar), ctx);
    return result;
  }

  case Expr::UDiv: 
  case Expr::SDiv: {
    SignHandler sh(e->getKind() == Expr::SDiv);
    BinaryExpr*de = cast<BinaryExpr>(e);
    Z3ASTHandle left = construct(de->left, width_out);
    assert(*width_out != 1 && "uncanonicalized sdiv");
    Z3ASTHandle right = construct(de->right, width_out);
    Z3ASTHandle result = Z3ASTHandle(Z3_mk_div(ctx, left, right), ctx);
    return result;
  }

  case Expr::URem: 
  case Expr::SRem: {
    SignHandler sh(e->getKind() == Expr::SRem);
    BinaryExpr*de = cast<BinaryExpr>(e);
    Z3ASTHandle left = construct(de->left, width_out);
    Z3ASTHandle right = construct(de->right, width_out);
    assert(*width_out != 1 && "uncanonicalized srem");
    // LLVM's srem instruction says that the sign follows the dividend
    // (``left``).
    // Z3's C API says ``Z3_mk_bvsrem()`` does this so these seem to match.
    Z3ASTHandle result = Z3ASTHandle(Z3_mk_rem(ctx, left, right), ctx);
    return result;
  }

  // Bitwise
  case Expr::Not: {
    NotExpr *ne = cast<NotExpr>(e);
    Z3ASTHandle expr = construct(ne->expr, width_out);
    if (*width_out == 1) {
      return notExpr(expr);
    } else {
      assert(0 && "Can't do bitwise not");
    }
  }

  case Expr::And: {
    AndExpr *ae = cast<AndExpr>(e);
    Z3ASTHandle left = construct(ae->left, width_out);
    Z3ASTHandle right = construct(ae->right, width_out);
    if (*width_out == 1) {
      return andExpr(left, right);
    } else {
      assert(0 && "Can't do bitwise and");
    }
  }

  case Expr::Or: {
    OrExpr *oe = cast<OrExpr>(e);
    Z3ASTHandle left = construct(oe->left, width_out);
    Z3ASTHandle right = construct(oe->right, width_out);
    if (*width_out == 1) {
      return orExpr(left, right);
    } else {
      assert(0 && "Can't do bitwise or");
    }
  }

  case Expr::Xor: {
    XorExpr *xe = cast<XorExpr>(e);
    Z3ASTHandle left = construct(xe->left, width_out);
    Z3ASTHandle right = construct(xe->right, width_out);

    if (*width_out == 1) {
      // XXX check for most efficient?
      return  Z3ASTHandle(Z3_mk_xor(ctx, left, right), ctx);
    } else {
      assert(0 && "Can't do bitwise xor");
    }
  }

  case Expr::Shl: 
  case Expr::LShr: 
  case Expr::AShr: {
      assert(0 && "Int solver doesn't do shifts");
  }

  // Comparison

  case Expr::Eq: {
    EqExpr *ee = cast<EqExpr>(e);
    Z3ASTHandle left = construct(ee->left, width_out);
    Z3ASTHandle right = construct(ee->right, width_out);
    if (*width_out == 1) {
      if (ConstantExpr *CE = dyn_cast<ConstantExpr>(ee->left)) {
        if (CE->isTrue())
          return right;
        return notExpr(right);
      } else {
        return iffExpr(left, right);
      }
    } else {
      *width_out = 1;
      return eqExpr(left, right);
    }
  }

  case Expr::Ult: 
  case Expr::Slt: {
    SignHandler sh(e->getKind() == Expr::Slt);
    BinaryExpr *se = cast<BinaryExpr>(e);
    Z3ASTHandle left = construct(se->left, width_out);
    Z3ASTHandle right = construct(se->right, width_out);
    assert(*width_out != 1 && "uncanonicalized slt");
    *width_out = 1;
    return Z3ASTHandle(Z3_mk_lt(ctx, left, right), ctx);
  }

  case Expr::Ule: 
  case Expr::Sle: {
    SignHandler sh(e->getKind() == Expr::Sle);
    BinaryExpr*se = cast<BinaryExpr>(e);
    Z3ASTHandle left = construct(se->left, width_out);
    Z3ASTHandle right = construct(se->right, width_out);
    assert(*width_out != 1 && "uncanonicalized sle");
    *width_out = 1;
    return Z3ASTHandle(Z3_mk_le(ctx, left, right), ctx);
  }

// unused due to canonicalization
#if 0
  case Expr::Ne:
  case Expr::Ugt:
  case Expr::Uge:
  case Expr::Sgt:
  case Expr::Sge:
#endif

  default:
    assert(0 && "unhandled Expr type");
    return getTrue();
  }
}
}
#endif // ENABLE_Z3
