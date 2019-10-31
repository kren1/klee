//===-- ArrayExprOptimizer.cpp --------------------------------------------===//
//
//                     The KLEE Symbolic Virtual Machine
//
// This file is distributed under the University of Illinois Open Source
// License. See LICENSE.TXT for details.
//
//===----------------------------------------------------------------------===//

#include "ArrayExprOptimizer.h"
#include "ArrayExprRewriter.h"
#include "ArrayExprVisitor.h"
#include "AssignmentGenerator.h"

#include "klee/Config/Version.h"
#include "klee/Expr/Assignment.h"
#include "klee/Expr/ArrayRanges.h"
#include "klee/Expr/ExprBuilder.h"
#include "klee/Internal/Support/ErrorHandling.h"
#include "klee/OptionCategories.h"
#include "klee/util/BitArray.h"

#include <llvm/ADT/APInt.h>
#include <llvm/Support/Casting.h>
#include <llvm/Support/CommandLine.h>

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <set>
#include <unordered_map>

using namespace klee;

namespace klee {
llvm::cl::opt<ArrayOptimizationType> OptimizeArray(
    "optimize-array",
    llvm::cl::values(clEnumValN(ALL, "all",
                                "Combining index and value transformations"),
                     clEnumValN(INDEX, "index", "Index-based transformation"),
                     clEnumValN(VALUE, "value",
                                "Value-based transformation at branch (both "
                                "concrete and concrete/symbolic)")
                         KLEE_LLVM_CL_VAL_END),
    llvm::cl::init(NONE),
    llvm::cl::desc("Optimize accesses to either concrete or concrete/symbolic "
                   "arrays. (default=false)"),
    llvm::cl::cat(klee::SolvingCat));

llvm::cl::opt<double> ArrayValueRatio(
    "array-value-ratio",
    llvm::cl::desc("Maximum ratio of unique values to array size for which the "
                   "value-based transformations are applied."),
    llvm::cl::init(1.0), llvm::cl::value_desc("Unique Values / Array Size"),
    llvm::cl::cat(klee::SolvingCat));

llvm::cl::opt<double> ArrayValueSymbRatio(
    "array-value-symb-ratio",
    llvm::cl::desc("Maximum ratio of symbolic values to array size for which "
                   "the mixed value-based transformations are applied."),
    llvm::cl::init(1.0), llvm::cl::value_desc("Symbolic Values / Array Size"),
    llvm::cl::cat(klee::SolvingCat));
}; // namespace klee

ref<Expr> extendRead(const UpdateList &ul, const ref<Expr> index,
                     Expr::Width w) {
  switch (w) {
  default:
    assert(0 && "invalid width");
  case Expr::Int8:
    return ReadExpr::alloc(ul, index);
  case Expr::Int16:
    return ConcatExpr::create(
        ReadExpr::alloc(
            ul, AddExpr::create(ConstantExpr::create(1, Expr::Int32), index)),
        ReadExpr::alloc(ul, index));
  case Expr::Int32:
    return ConcatExpr::create4(
        ReadExpr::alloc(
            ul, AddExpr::create(ConstantExpr::create(3, Expr::Int32), index)),
        ReadExpr::alloc(
            ul, AddExpr::create(ConstantExpr::create(2, Expr::Int32), index)),
        ReadExpr::alloc(
            ul, AddExpr::create(ConstantExpr::create(1, Expr::Int32), index)),
        ReadExpr::alloc(ul, index));
  case Expr::Int64:
    return ConcatExpr::create8(
        ReadExpr::alloc(
            ul, AddExpr::create(ConstantExpr::create(7, Expr::Int32), index)),
        ReadExpr::alloc(
            ul, AddExpr::create(ConstantExpr::create(6, Expr::Int32), index)),
        ReadExpr::alloc(
            ul, AddExpr::create(ConstantExpr::create(5, Expr::Int32), index)),
        ReadExpr::alloc(
            ul, AddExpr::create(ConstantExpr::create(4, Expr::Int32), index)),
        ReadExpr::alloc(
            ul, AddExpr::create(ConstantExpr::create(3, Expr::Int32), index)),
        ReadExpr::alloc(
            ul, AddExpr::create(ConstantExpr::create(2, Expr::Int32), index)),
        ReadExpr::alloc(
            ul, AddExpr::create(ConstantExpr::create(1, Expr::Int32), index)),
        ReadExpr::alloc(ul, index));
  }
}

ref<Expr> ExprOptimizer::optimizeExpr(const ref<Expr> &e, bool valueOnly) {
  // Nothing to optimise for constant expressions
  if (isa<ConstantExpr>(e))
    return e;

  // If no is optimization enabled, return early avoid cache lookup
  if (OptimizeArray == NONE)
    return e;

  unsigned hash = e->hash();
  if (cacheExprUnapplicable.find(hash) != cacheExprUnapplicable.end())
    return e;

  // Find cached expressions
  auto cached = cacheExprOptimized.find(hash);
  if (cached != cacheExprOptimized.end())
    return cached->second;

  ref<Expr> result;
  // ----------------------- INDEX-BASED OPTIMIZATION -------------------------
  if (!valueOnly && (OptimizeArray == ALL || OptimizeArray == INDEX)) {
    array2idx_ty arrays;
    ConstantArrayExprVisitor aev(arrays);
    aev.visit(e);

    if (arrays.empty() || aev.isIncompatible()) {
      // We do not optimize expressions other than those with concrete
      // arrays with a symbolic index
      // If we cannot optimize the expression, we return a failure only
      // when we are not combining the optimizations
      if (OptimizeArray == INDEX) {
        cacheExprUnapplicable.insert(hash);
        return e;
      }
    } else {
      mapIndexOptimizedExpr_ty idx_valIdx;

      // Compute those indexes s.t. orig_expr =equisat= (k==i|k==j|..)
      if (computeIndexes(arrays, e, idx_valIdx)) {
        if (!idx_valIdx.empty()) {
          // Create new expression on indexes
          result = ExprRewriter::createOptExpr(e, arrays, idx_valIdx);
        } else {
          klee_warning("OPT_I: infeasible branch!");
          result = ConstantExpr::create(0, Expr::Bool);
        }
        // Add new expression to cache
        if (result.get()) {
          klee_warning("OPT_I: successful");
          cacheExprOptimized[hash] = result;
        } else {
          klee_warning("OPT_I: unsuccessful");
        }
      } else {
        klee_warning("OPT_I: unsuccessful");
        cacheExprUnapplicable.insert(hash);
      }
    }
  }
  // ----------------------- VALUE-BASED OPTIMIZATION -------------------------
  if (OptimizeArray == VALUE ||
      (OptimizeArray == ALL && (!result.get() || valueOnly))) {
    std::vector<const ReadExpr *> reads;
    std::map<const ReadExpr *, std::pair<unsigned, Expr::Width>> readInfo;
    ArrayReadExprVisitor are(reads, readInfo);
    are.visit(e);
    std::reverse(reads.begin(), reads.end());

    if (reads.empty() || are.isIncompatible()) {
      cacheExprUnapplicable.insert(hash);
      return e;
    }

    ref<Expr> selectOpt =
        getSelectOptExpr(e, reads, readInfo, are.containsSymbolic());
    if (selectOpt.get()) {
      klee_warning("OPT_V: successful");
      result = selectOpt;
      cacheExprOptimized[hash] = result;
    } else {
      klee_warning("OPT_V: unsuccessful");
      cacheExprUnapplicable.insert(hash);
    }
  }
  if (result.isNull())
    return e;
  return result;
}

bool ExprOptimizer::computeIndexes(array2idx_ty &arrays, const ref<Expr> &e,
                                   mapIndexOptimizedExpr_ty &idx_valIdx) const {
  bool success = false;
  // For each constant array found
  for (auto &element : arrays) {
    const Array *arr = element.first;

    assert(arr->isConstantArray() && "Array is not concrete");
    assert(element.second.size() == 1 && "Multiple indexes on the same array");

    IndexTransformationExprVisitor idxt_v(arr);
    idxt_v.visit(e);
    assert((idxt_v.getWidth() % arr->range == 0) && "Read is not aligned");
    Expr::Width width = idxt_v.getWidth() / arr->range;

    if (idxt_v.getMul().get()) {
      // If we have a MulExpr in the index, we can optimize our search by
      // skipping all those indexes that are not multiple of such value.
      // In fact, they will be rejected by the MulExpr interpreter since it
      // will not find any integer solution
      auto e = idxt_v.getMul();
      auto ce = dyn_cast<ConstantExpr>(e);
      assert(ce && "Not a constant expression");
      uint64_t mulVal = (*ce->getAPValue().getRawData());
      // So far we try to limit this optimization, but we may try some more
      // aggressive conditions (i.e. mulVal > width)
      if (width == 1 && mulVal > 1)
        width = mulVal;
    }

    // For each concrete value 'i' stored in the array
    for (size_t aIdx = 0; aIdx < arr->constantValues.size(); aIdx += width) {
      auto *a = new Assignment();
      std::vector<const Array *> objects;
      std::vector<std::vector<unsigned char>> values;

      // For each symbolic index Expr(k) found
      for (auto &index_it : element.second) {
        ref<Expr> idx = index_it;
        ref<Expr> val = ConstantExpr::alloc(aIdx, arr->getDomain());
        // We create a partial assignment on 'k' s.t. Expr(k)==i
        bool assignmentSuccess =
            AssignmentGenerator::generatePartialAssignment(idx, val, a);
        success |= assignmentSuccess;

        // If the assignment satisfies both the expression 'e' and the PC
        ref<Expr> evaluation = a->evaluate(e);
        if (assignmentSuccess && evaluation->isTrue()) {
          if (idx_valIdx.find(idx) == idx_valIdx.end()) {
            idx_valIdx.insert(std::make_pair(idx, std::vector<ref<Expr>>()));
          }
          idx_valIdx[idx].emplace_back(
              ConstantExpr::alloc(aIdx, arr->getDomain()));
        }
      }
      delete a;
    }
  }
  return success;
}

ref<Expr> ExprOptimizer::getSelectOptExpr(
    const ref<Expr> &e, std::vector<const ReadExpr *> &reads,
    std::map<const ReadExpr *, std::pair<unsigned, Expr::Width>> &readInfo,
    bool isSymbolic) {
  ref<Expr> notFound;
  ref<Expr> toReturn;

  // Array is concrete
  if (!isSymbolic) {
    std::map<unsigned, ref<Expr>> optimized;
    for (auto &read : reads) {
      auto info = readInfo[read];
      auto cached = cacheReadExprOptimized.find(read->hash());
      if (cached != cacheReadExprOptimized.end()) {
        optimized.insert(std::make_pair(info.first, (*cached).second));
        continue;
      }
      Expr::Width width = read->getWidth();
      if (info.second > width) {
        width = info.second;
      }
      unsigned size = read->updates.root->getSize();
      unsigned bytesPerElement = width / 8;
      unsigned elementsInArray = size / bytesPerElement;

      // Note: we already filtered the ReadExpr, so here we can safely
      // assume that the UpdateNodes contain ConstantExpr indexes and values
      assert(read->updates.root->isConstantArray() &&
             "Expected concrete array, found symbolic array");
      auto arrayConstValues = read->updates.root->constantValues;
      for (const UpdateNode *un = read->updates.head; un; un = un->next) {
        auto ce = dyn_cast<ConstantExpr>(un->index);
        assert(ce && "Not a constant expression");
        uint64_t index = ce->getAPValue().getZExtValue();
        assert(index < arrayConstValues.size());
        auto arrayValue = dyn_cast<ConstantExpr>(un->value);
        assert(arrayValue && "Not a constant expression");
        arrayConstValues[index] = arrayValue;
      }
      std::vector<uint64_t> arrayValues;
      // Get the concrete values from the array
      for (unsigned i = 0; i < elementsInArray; i++) {
        uint64_t val = 0;
        for (unsigned j = 0; j < bytesPerElement; j++) {
          val |= (*(
                       arrayConstValues[(i * bytesPerElement) + j]
                           .get()
                           ->getAPValue()
                           .getRawData())
                  << (j * 8));
        }
        arrayValues.push_back(val);
      }

      ref<Expr> index = read->index;
      IndexCleanerVisitor ice;
      index = ice.visit(index);
     // index->dump();
////      llvm::errs() << read->updates.root->getName() << ": ";
//      llvm::errs() << "width: " << width << "\n";
//         const Array* A = read->updates.root;
//       llvm::errs() << A->getName() << ": [";
//        for (unsigned i = 0, e = A->size; i != e; ++i) {
//          if (i)
//            llvm::errs() << " ";
//          llvm::errs() << A->constantValues[i];
//        }
//        llvm::errs() << "]\n";


      ref<Expr> opt =
          buildConstantSelectExpr(index, arrayValues, width, elementsInArray);
      if (opt.get()) {
        cacheReadExprOptimized[read->hash()] = opt;
        optimized.insert(std::make_pair(info.first, opt));
      }
    }
    ArrayValueOptReplaceVisitor replacer(optimized);
    toReturn = replacer.visit(e);
  }

  // Array is mixed concrete/symbolic
  // \pre: array is concrete && updatelist contains at least one symbolic value
  // OR
  //       array is symbolic && updatelist contains at least one concrete value
  else {
    std::map<unsigned, ref<Expr>> optimized;
    for (auto &read : reads) {
      auto info = readInfo[read];
      auto cached = cacheReadExprOptimized.find(read->hash());
      if (cached != cacheReadExprOptimized.end()) {
        optimized.insert(std::make_pair(info.first, (*cached).second));
        continue;
      }
      Expr::Width width = read->getWidth();
      if (info.second > width) {
        width = info.second;
      }
      unsigned size = read->updates.root->getSize();
      unsigned bytesPerElement = width / 8;
      unsigned elementsInArray = size / bytesPerElement;
      bool symbArray = read->updates.root->isSymbolicArray();

      BitArray ba(size, symbArray);
      // Note: we already filtered the ReadExpr, so here we can safely
      // assume that the UpdateNodes contain ConstantExpr indexes, but in
      // this case we *cannot* assume anything on the values
      auto arrayConstValues = read->updates.root->constantValues;
      if (arrayConstValues.size() < size) {
        // We need to "force" initialization of the values
        for (size_t i = arrayConstValues.size(); i < size; i++) {
          arrayConstValues.push_back(ConstantExpr::create(0, Expr::Int8));
        }
      }
      for (const UpdateNode *un = read->updates.head; un; un = un->next) {
        auto ce = dyn_cast<ConstantExpr>(un->index);
        assert(ce && "Not a constant expression");
        uint64_t index = ce->getAPValue().getLimitedValue();
        if (!isa<ConstantExpr>(un->value)) {
          ba.set(index);
        } else {
          ba.unset(index);
          auto arrayValue =
              dyn_cast<ConstantExpr>(un->value);
          assert(arrayValue && "Not a constant expression");
          arrayConstValues[index] = arrayValue;
        }
      }

      std::vector<std::pair<uint64_t, bool>> arrayValues;
      unsigned symByteNum = 0;
      for (unsigned i = 0; i < elementsInArray; i++) {
        uint64_t val = 0;
        bool elementIsConcrete = true;
        for (unsigned j = 0; j < bytesPerElement; j++) {
          if (ba.get((i * bytesPerElement) + j)) {
            elementIsConcrete = false;
            break;
          } else {
            val |= (*(
                         arrayConstValues[(i * bytesPerElement) + j]
                             .get()
                             ->getAPValue()
                             .getRawData())
                    << (j * 8));
          }
        }
        if (elementIsConcrete) {
          arrayValues.emplace_back(val, true);
        } else {
          symByteNum++;
          arrayValues.emplace_back(0, false);
        }
      }

      if (((double)symByteNum / (double)elementsInArray) <=
          ArrayValueSymbRatio) {
        // If the optimization can be applied we apply it
        // Build the dynamic select expression
        ref<Expr> opt =
            buildMixedSelectExpr(read, arrayValues, width, elementsInArray);
        if (opt.get()) {
          cacheReadExprOptimized[read->hash()] = opt;
          optimized.insert(std::make_pair(info.first, opt));
        }
      }
    }
    ArrayValueOptReplaceVisitor replacer(optimized, false);
    toReturn = replacer.visit(e);
  }

  return toReturn.get() ? toReturn : notFound;
}

ref<Expr> ExprOptimizer::buildConstantSelectExpr(
    const ref<Expr> &index, std::vector<uint64_t> &arrayValues,
    Expr::Width width, unsigned arraySize) const {
//        llvm::errs() << "Constant select expr\n";
  std::set<uint64_t> unique_array_values;
  ExprBuilder *builder =createDefaultExprBuilder();
  builder = createConstantFoldingExprBuilder(builder);
  Expr::Width valWidth = width;
  ref<Expr> result;

  ref<Expr> actualIndex;
  if (index->getWidth() > Expr::Int32) {
    actualIndex = ExtractExpr::alloc(index, 0, Expr::Int32);
  } else {
    actualIndex = index;
  }
  Expr::Width idxWidth = actualIndex->getWidth();

  auto derivative = ArrayRanges::firstDerivative(arrayValues);
  auto derivativeRanges = ArrayRanges::equalRanges(derivative);
  auto ranges = ArrayRanges::equalRanges(arrayValues);
 // llvm::errs() << "ranges: " << ranges << "\nderivatives[";
 // for(const auto& d : derivative)
 //     llvm::errs() << d << ", ";
 // llvm::errs() << "]\n";

  if (((double)ranges.size() / (double)(arraySize + 1)) >=
      ArrayValueRatio) {
    return result;
  }

  std::unordered_map<uint64_t, std::vector<ArrayRanges::Range>> exprMap;
  for (const auto& range : ranges) {
    auto value = arrayValues[range.start];
    if (exprMap.count(value) > 0) {
      exprMap[value].emplace_back(range);
    } else {
      exprMap[value] = {range};
    }
  }
  for (const auto& valRange : exprMap) {
  //    llvm::errs() << valRange.first << " -> " << valRange.second << "\n";
  }

 //   llvm::errs() << "derivative Ranges: " << derivativeRanges << "\n";
  std::vector<std::pair<LinFun, ArrayRanges::Range>> derExprMap;
  assert(derivativeRanges.size() > 1 && "TODO edge case");

  auto rangeIt = derivativeRanges.begin();
  assert(rangeIt->isUnitary() && "First value should always be a transition (ie. infinite derivative)"); 
  auto arrayVal = arrayValues[rangeIt->start];

  auto fun = LinFun(arrayVal,1, rangeIt->start);
  assert(fun.eval(rangeIt->start) == arrayVal);
  derExprMap.emplace_back(fun, *rangeIt);
  rangeIt++;
//  llvm::errs() << derExprMap.back().second << " -> " << derExprMap.back().first << "\n";


  while(rangeIt != derivativeRanges.end()) {
    arrayVal = arrayValues[rangeIt->start];
    if(
        //(derivative[rangeIt->start] != derivative[derExprMap.back().second.start])
       (derExprMap.back().second.isUnitary() && !rangeIt->isUnitary())) { //linear case
      auto& record = derExprMap.back();
      record.first.coef = derivative[rangeIt->start];
      assert(record.second.merge(*rangeIt));
    } else { //transition case
      fun = LinFun(arrayVal, 1, rangeIt->start);
      derExprMap.emplace_back(fun, *rangeIt);
    }
//    llvm::errs() << derExprMap.back().second << " -> " << derExprMap.back().first << "\n";
//    llvm::errs() << arrayVal << " x: " << rangeIt->start << "\n";
    assert(derExprMap.back().first.eval(rangeIt->start) == arrayVal);
    rangeIt++;
  }

  int ct = 0;
#define valAPInt(value) llvm::APInt(valWidth, (uint64_t)value, false)
#define idxAPInt(value) llvm::APInt(idxWidth, (uint64_t)value, false)
#define idxsAPInt(value) llvm::APInt(idxWidth, (int64_t)value, true)
  auto createIdxForRange = [&](const ArrayRanges::Range& range) {
        return range.isUnitary() 
              ?
                EqExpr::create(actualIndex,
                             builder->Constant(idxAPInt(range.start)))

              : AndExpr::create(
                  SgeExpr::create(actualIndex,
                                  builder->Constant(idxAPInt(range.start))),
                  SltExpr::create(
                      actualIndex,
                      builder->Constant(idxAPInt(range.end))));
  };
  auto linFunToExpr = [&](const LinFun& fun, const ref<Expr>& idx) {
      ref<Expr> ret = builder->Sub(idx, builder->Constant(idxsAPInt(fun.x_offset)));
      ret = builder->Mul(ret, builder->Constant(idxsAPInt(fun.coef)));
      ret = builder->Add(ret, builder->Constant(idxsAPInt(fun.base)));
      return ret;
  };
  if(derExprMap.size() < exprMap.size()) {
      klee_warning("Using linear regions");
      for (const auto& linFunRange : derExprMap) {
        auto linFun = linFunRange.first;
        auto range = linFunRange.second;
  //      llvm::errs() << range<< " -> " << linFun << "\n";
        ref<Expr> temp;
        if (ct == 0) {
          temp = linFunToExpr(linFun, actualIndex);
        } else {
          temp = SelectExpr::create(createIdxForRange(range),
                      linFunToExpr(linFun, actualIndex),
                      result);
          
        }
        result = temp;
        ct++;
      }
  } else {
      klee_warning("Using constant regions");
      for (const auto& range : exprMap) {
        auto value = range.first;
        auto rangesForValue = range.second;
        ref<Expr> temp;
        if (ct == 0) {
          temp = builder->Constant(valAPInt(value));
        } else {
          if (rangesForValue.size() == 1) {
            temp = SelectExpr::create(createIdxForRange(rangesForValue[0]),
                      builder->Constant(valAPInt(value)),
                      result);

          } else {
            ref<Expr> currOr;
            currOr = createIdxForRange(rangesForValue[0]);

            for (size_t i = 1; i < rangesForValue.size(); i++) {
              ref<Expr> tempOr;
              tempOr = OrExpr::create( 
                createIdxForRange(rangesForValue[i]),
                currOr
              );
              currOr = tempOr;
            }
            temp = SelectExpr::create(currOr, builder->Constant(valAPInt(value)),
                                      result);
          }
        }
        result = temp;
        ct++;
      }
  }

  delete (builder);

  return result;
}

ref<Expr> ExprOptimizer::buildMixedSelectExpr(
    const ReadExpr *re, std::vector<std::pair<uint64_t, bool>> &arrayValues,
    Expr::Width width, unsigned elementsInArray) const {
  ExprBuilder *builder = createDefaultExprBuilder();
  std::vector<uint64_t> values;
  std::vector<std::pair<uint64_t, uint64_t>> ranges;
  std::vector<uint64_t> holes;
  std::set<uint64_t> unique_array_values;

  unsigned arraySize = elementsInArray;
  unsigned curr_idx = 0;
  uint64_t curr_val = arrayValues[0].first;

  bool emptyRange = true;
  // Calculate Range values
  for (size_t i = 0; i < arrayValues.size(); i++) {
    // If the value is concrete
    if (arrayValues[i].second) {
      // The range contains a concrete value
      emptyRange = false;
      uint64_t temp = arrayValues[i].first;
      unique_array_values.insert(temp);
      if (temp != curr_val) {
        ranges.emplace_back(curr_idx, i);
        values.push_back(curr_val);
        curr_val = temp;
        curr_idx = i;
        if (i == (arraySize - 1)) {
          ranges.emplace_back(curr_idx, curr_idx + 1);
          values.push_back(curr_val);
        }
      } else if (i == (arraySize - 1)) {
        ranges.emplace_back(curr_idx, i + 1);
        values.push_back(curr_val);
      }
    } else {
      holes.push_back(i);
      // If this is not an empty range
      if (!emptyRange) {
        ranges.emplace_back(curr_idx, i);
        values.push_back(curr_val);
      }
      curr_val = arrayValues[i + 1].first;
      curr_idx = i + 1;
      emptyRange = true;
    }
  }

  assert(!unique_array_values.empty() && "No unique values");
  assert(!ranges.empty() && "No ranges");

  ref<Expr> result;
  if (((double)unique_array_values.size() / (double)(arraySize)) <=
      ArrayValueRatio) {
    // The final "else" expression will be the original unoptimized array read
    // expression
    unsigned range_start = 0;
    if (holes.empty()) {
      result = builder->Constant(llvm::APInt(width, values[0], false));
      range_start = 1;
    } else {
      ref<Expr> firstIndex = MulExpr::create(
          ConstantExpr::create(holes[0], re->index->getWidth()),
          ConstantExpr::create(width / 8, re->index->getWidth()));
      result = extendRead(re->updates, firstIndex, width);
      for (size_t i = 1; i < holes.size(); i++) {
        ref<Expr> temp_idx = MulExpr::create(
            ConstantExpr::create(holes[i], re->index->getWidth()),
            ConstantExpr::create(width / 8, re->index->getWidth()));
        ref<Expr> cond = EqExpr::create(
            re->index, ConstantExpr::create(holes[i], re->index->getWidth()));
        ref<Expr> temp = SelectExpr::create(
            cond, extendRead(re->updates, temp_idx, width), result);
        result = temp;
      }
    }

    ref<Expr> new_index = re->index;
    IndexCleanerVisitor ice;
    new_index = ice.visit(new_index);

    int new_index_width = new_index->getWidth();
    // Iterate through all the ranges
    for (size_t i = range_start; i < ranges.size(); i++) {
      ref<Expr> temp;
      if (ranges[i].second - 1 == ranges[i].first) {
        ref<Expr> cond = EqExpr::create(
            new_index, ConstantExpr::create(ranges[i].first, new_index_width));
        ref<Expr> t = ConstantExpr::create(values[i], width);
        ref<Expr> f = result;
        temp = SelectExpr::create(cond, t, f);
      } else {
        // Create the select constraint
        ref<Expr> cond = AndExpr::create(
            SgeExpr::create(new_index, ConstantExpr::create(ranges[i].first,
                                                            new_index_width)),
            SltExpr::create(new_index, ConstantExpr::create(ranges[i].second,
                                                            new_index_width)));
        ref<Expr> t = ConstantExpr::create(values[i], width);
        ref<Expr> f = result;
        temp = SelectExpr::create(cond, t, f);
      }
      result = temp;
    }
  }

  delete (builder);

  return result;
}
