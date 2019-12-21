//===-- CachingSolver.cpp - Caching expression solver ---------------------===//
//
//                     The KLEE Symbolic Virtual Machine
//
// This file is distributed under the University of Illinois Open Source
// License. See LICENSE.TXT for details.
//
//===----------------------------------------------------------------------===//


#include "klee/Solver/Solver.h"

#include "klee/Expr/Constraints.h"
#include "klee/Expr/Expr.h"
#include "klee/util/WeakRef.h"
#include "klee/Solver/IncompleteSolver.h"
#include "klee/Solver/SolverImpl.h"
#include "klee/Solver/SolverStats.h"

#include "llvm/Support/CommandLine.h"

#include <unordered_map>
#define unordered_map std::unordered_map
//#include <ciso646>
//#ifdef _LIBCPP_VERSION
//#define unordered_map std::unordered_map
//#else
//#include <tr1/unordered_map>
//#define unordered_map std::tr1::unordered_map
//#endif

using namespace klee;
  llvm::cl::opt<unsigned> CacheSize(
    "cache-size", llvm::cl::init(100),
    llvm::cl::desc("Size of cache (default=100)"));



SQLIntStatistic oneRefEntries("CacheOneRefEntries", "COneR");
SQLIntStatistic oneConstraintOneRefEntries("CacheOneConstraintRefEntries", "COneCR");
SQLIntStatistic cacheSize("BranchCacheSize", "BCachS");
class CachingSolver : public SolverImpl {
private:
  ref<Expr> canonicalizeQuery(ref<Expr> originalQuery,
                              bool &negationUsed);

  void cacheInsert(const Query& query,
                   IncompleteSolver::PartialValidity result);

  bool cacheLookup(const Query& query,
                   IncompleteSolver::PartialValidity &result);
  
  struct CacheEntry {
    std::vector<weak_ref<Expr>> key;
    CacheEntry(const ConstraintManager &c, ref<Expr> &q) {
        assert(q.get());
        key.emplace_back(q.get());
        for(const auto& con: c) {
            assert(con.get());
            key.emplace_back(con.get());
        }
    }


    CacheEntry(const CacheEntry &ce)
      : key(ce.key) {}
    

    bool operator==(const CacheEntry &b) const {
      bool isNull = false;
      for(auto& e : key) {
          if(e.isNull()) {
              isNull = true;
              break;
          }
      }
      if(isNull) {
          for(auto& e : b.key) {
              if(e.isNull()) return true;
          }
      }
      if(isNull) return false;
      return key == b.key;
    }
  };

  struct CacheEntryHash {
    unsigned operator()(const CacheEntry &ce) const {
      unsigned result = 0;

      for (auto const &constraint : ce.key) {
        if(constraint.isNull()) return 0;
        result ^= constraint->hash();
      }

      return result;
    }
  };

  typedef unordered_map<CacheEntry, 
                        IncompleteSolver::PartialValidity, 
                        CacheEntryHash> cache_map;
  
  Solver *solver;
  cache_map cache;

public:
  CachingSolver(Solver *s) : solver(s) { cache.reserve(CacheSize);}
  ~CachingSolver() { cache.clear(); delete solver; }

  static bool check;
  void countRefCounts() {

      static int prevOneRef = 0;
      static int prevSize = 0;
      static int prevoneConstraintOneRef = 0;
      if(!check) return;
      int oneRefs = 0;
      int oneConstraintOneRef = 0;
      for(auto& e : cache) {
          auto& entry = e.first;
          assert(entry.key.size() > 0);
          bool allZeroRefs =  true;
          bool oneZeroRef = false;
          for(const auto& constraint : entry.key) {

              int cnt = constraint.isNull() ? 0 : constraint->_refCount.refCount;
              if(cnt == 0) { 
//                  llvm::errs() << "\tZero strong refs, weak refs: " << constraint->_refCount.weakrefCount << "\n";
              }
              allZeroRefs &= cnt == 0;
              oneZeroRef |= cnt == 0;
          }
          if(allZeroRefs) oneRefs++;
          if(oneZeroRef) oneConstraintOneRef++;

      }
      oneRefEntries += oneRefs - prevOneRef;
      oneConstraintOneRefEntries += oneConstraintOneRef - prevoneConstraintOneRef;
      cacheSize += cache.size() - prevSize;
      llvm::errs() << "============ One refs: " << oneRefs << " oneConstraintZeroRef: " << oneConstraintOneRef <<  " cache size: " << cache.size() <<  "\n";
      llvm::errs() << "load factor: " << cache.load_factor() << " max_load factor: " << cache.max_load_factor() << " max size: " << cache.max_size() <<  "\n";
      prevOneRef = oneRefs;
      prevoneConstraintOneRef = oneConstraintOneRef;
      prevSize = cache.size();
      check = false;
  }

  bool computeValidity(const Query&, Solver::Validity &result);
  bool computeTruth(const Query&, bool &isValid);
  bool computeValue(const Query& query, ref<Expr> &result) {
    ++stats::queryCacheMisses;
    return solver->impl->computeValue(query, result);
  }
  bool computeInitialValues(const Query& query,
                            const std::vector<const Array*> &objects,
                            std::vector< std::vector<unsigned char> > &values,
                            bool &hasSolution) {
    ++stats::queryCacheMisses;
    return solver->impl->computeInitialValues(query, objects, values, 
                                              hasSolution);
  }
  SolverRunStatus getOperationStatusCode();
  char *getConstraintLog(const Query&);
  void setCoreSolverTimeout(time::Span timeout);
};

bool CachingSolver::check = false;
/** @returns the canonical version of the given query.  The reference
    negationUsed is set to true if the original query was negated in
    the canonicalization process. */
ref<Expr> CachingSolver::canonicalizeQuery(ref<Expr> originalQuery,
                                           bool &negationUsed) {
  ref<Expr> negatedQuery = Expr::createIsZero(originalQuery);

  // select the "smaller" query to the be canonical representation
  if (originalQuery.compare(negatedQuery) < 0) {
    negationUsed = false;
    return originalQuery;
  } else {
    negationUsed = true;
    return negatedQuery;
  }
}

/** @returns true on a cache hit, false of a cache miss.  Reference
    value result only valid on a cache hit. */
bool CachingSolver::cacheLookup(const Query& query,
                                IncompleteSolver::PartialValidity &result) {
  bool negationUsed;
  ref<Expr> canonicalQuery = canonicalizeQuery(query.expr, negationUsed);
  countRefCounts();

  CacheEntry ce(query.constraints, canonicalQuery);
  cache_map::iterator it = cache.find(ce);
  
  if (it != cache.end()) {
    result = (negationUsed ?
              IncompleteSolver::negatePartialValidity(it->second) :
              it->second);
    return true;
  }
  
  return false;
}

/// Inserts the given query, result pair into the cache.
void CachingSolver::cacheInsert(const Query& query,
                                IncompleteSolver::PartialValidity result) {
  if(cache.load_factor() > 0.95*cache.max_load_factor()) {
//      llvm::errs() << "Rehashing!\n";
      cache = cache_map(cache.begin(), cache.end(), cache.bucket_count()*2);
  }
  bool negationUsed;
  ref<Expr> canonicalQuery = canonicalizeQuery(query.expr, negationUsed);

  CacheEntry ce(query.constraints, canonicalQuery);
  IncompleteSolver::PartialValidity cachedResult = 
    (negationUsed ? IncompleteSolver::negatePartialValidity(result) : result);
  
  cache.insert(std::make_pair(ce, cachedResult));
}

bool CachingSolver::computeValidity(const Query& query,
                                    Solver::Validity &result) {
  IncompleteSolver::PartialValidity cachedResult;
  bool tmp, cacheHit = cacheLookup(query, cachedResult);
  
  if (cacheHit) {
    switch(cachedResult) {
    case IncompleteSolver::MustBeTrue:   
      result = Solver::True;
      ++stats::queryCacheHits;
      return true;
    case IncompleteSolver::MustBeFalse:  
      result = Solver::False;
      ++stats::queryCacheHits;
      return true;
    case IncompleteSolver::TrueOrFalse:  
      result = Solver::Unknown;
      ++stats::queryCacheHits;
      return true;
    case IncompleteSolver::MayBeTrue: {
      ++stats::queryCacheMisses;
      if (!solver->impl->computeTruth(query, tmp))
        return false;
      if (tmp) {
        cacheInsert(query, IncompleteSolver::MustBeTrue);
        result = Solver::True;
        return true;
      } else {
        cacheInsert(query, IncompleteSolver::TrueOrFalse);
        result = Solver::Unknown;
        return true;
      }
    }
    case IncompleteSolver::MayBeFalse: {
      ++stats::queryCacheMisses;
      if (!solver->impl->computeTruth(query.negateExpr(), tmp))
        return false;
      if (tmp) {
        cacheInsert(query, IncompleteSolver::MustBeFalse);
        result = Solver::False;
        return true;
      } else {
        cacheInsert(query, IncompleteSolver::TrueOrFalse);
        result = Solver::Unknown;
        return true;
      }
    }
    default: assert(0 && "unreachable");
    }
  }

  ++stats::queryCacheMisses;
  
  if (!solver->impl->computeValidity(query, result))
    return false;

  switch (result) {
  case Solver::True: 
    cachedResult = IncompleteSolver::MustBeTrue; break;
  case Solver::False: 
    cachedResult = IncompleteSolver::MustBeFalse; break;
  default: 
    cachedResult = IncompleteSolver::TrueOrFalse; break;
  }
  
  cacheInsert(query, cachedResult);
  return true;
}

bool CachingSolver::computeTruth(const Query& query,
                                 bool &isValid) {
  IncompleteSolver::PartialValidity cachedResult;
  bool cacheHit = cacheLookup(query, cachedResult);

  // a cached result of MayBeTrue forces us to check whether
  // a False assignment exists.
  if (cacheHit && cachedResult != IncompleteSolver::MayBeTrue) {
    ++stats::queryCacheHits;
    isValid = (cachedResult == IncompleteSolver::MustBeTrue);
    return true;
  }

  ++stats::queryCacheMisses;
  
  // cache miss: query solver
  if (!solver->impl->computeTruth(query, isValid))
    return false;

  if (isValid) {
    cachedResult = IncompleteSolver::MustBeTrue;
  } else if (cacheHit) {
    // We know a true assignment exists, and query isn't valid, so
    // must be TrueOrFalse.
    assert(cachedResult == IncompleteSolver::MayBeTrue);
    cachedResult = IncompleteSolver::TrueOrFalse;
  } else {
    cachedResult = IncompleteSolver::MayBeFalse;
  }
  
  cacheInsert(query, cachedResult);
  return true;
}

SolverImpl::SolverRunStatus CachingSolver::getOperationStatusCode() {
  return solver->impl->getOperationStatusCode();
}

char *CachingSolver::getConstraintLog(const Query& query) {
  return solver->impl->getConstraintLog(query);
}

void CachingSolver::setCoreSolverTimeout(time::Span timeout) {
  solver->impl->setCoreSolverTimeout(timeout);
}

///

Solver *klee::createCachingSolver(Solver *_solver) {
  return new Solver(new CachingSolver(_solver));
}
