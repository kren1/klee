//===-- ArrayCache.h --------------------------------------------*- C++ -*-===//
//
//                     The KLEE Symbolic Virtual Machine
//
// This file is distributed under the University of Illinois Open Source
// License. See LICENSE.TXT for details.
//
//===----------------------------------------------------------------------===//

#ifndef KLEE_ARRAYRANGE_H
#define KLEE_ARRAYRANGE_H

#include <vector>
#include "klee/Expr/Expr.h"


namespace klee {

class ArrayRanges {
public:
    struct Range {
        uint32_t start; //inclusive
        uint32_t end; //not inclusive

        Range(uint32_t s, uint32_t e): start(s), end(e) {}

        //Does this range represent only 1 value
        inline bool isUnitary() const { return end - start == 1; }

        inline bool operator==(const Range& lhs) const {return start == lhs.start && end == lhs.end;}
        
        //Merges other into this range, returns true if successful
        inline bool merge(const Range& other) {
            if(end != other.start) return false;
            end = other.end;
            return true;
        }
    };
  template <class Arr> static std::vector<Range> equalRanges(const Arr &arr) {
    if(arr.size() == 0) return {};
    auto it = arr.begin();
    int rangeStart = 0;
    auto rangeValue = *it;
    int idx = 0;
    std::vector<Range> res;
    while(it != arr.end()) {
//        llvm::errs() << idx << ": " << rangeValue << " it:" << *it << "\n";
        if(rangeValue != *it) {
//            llvm::errs() << "emplace!\n";
            res.emplace_back(rangeStart, idx);
            rangeValue = *it;
            rangeStart = idx;
        }
        it++;
        idx++;
    }
    res.emplace_back(rangeStart, idx);
    return res;
  }
  template <class Arr> static std::vector<int> firstDerivative(const Arr &arr) {
    if(arr.size() < 1) assert(false && "TODO small derivatives");
    std::vector<int> res;
    res.reserve(arr.size());
    int prev =  2147483647; //Something big, so derivative is large
    for(const auto& elem : arr) {
        res.emplace_back(elem - prev);
        prev = elem;
    }
    return res;
  }
};

struct LinFun {
      int32_t base;
      int32_t coef;
      int32_t x_offset;
      LinFun(int32_t b, int32_t c, int32_t o): base(b), coef(c), x_offset(o) {}
      LinFun(int32_t b, int32_t c): base(b), coef(c), x_offset(0) {}
      LinFun(int32_t b): base(b), coef(0), x_offset(0) {}
      inline uint64_t eval(uint64_t val) const {return coef*(val - x_offset) + base; }
};


llvm::raw_ostream &operator<<(llvm::raw_ostream &os, const LinFun& r) {
//    if(r.coef == 0) return os << r.base;
//    if(r.base == 0) return os << r.coef << " * val";
//    if(r.coef == 1) return os << "val + " << r.base;
    return os <<  r.coef << " * (x - " << r.x_offset << ") + " << r.base;
}
llvm::raw_ostream &operator<<(llvm::raw_ostream &os, const ArrayRanges::Range& r) {
    return os << "[" << r.start << ", " << r.end << ")";
}
llvm::raw_ostream &operator<<(llvm::raw_ostream &os, const std::vector<ArrayRanges::Range>& vec) {
    os << "[";
    for(const auto& r : vec) {
        os << r << ", ";
    }
    os << "]";
    return os;
}
}

#undef unordered_set

#endif /* KLEE_ARRAYCACHE_H */
