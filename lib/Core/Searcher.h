//===-- Searcher.h ----------------------------------------------*- C++ -*-===//
//
//                     The KLEE Symbolic Virtual Machine
//
// This file is distributed under the University of Illinois Open Source
// License. See LICENSE.TXT for details.
//
//===----------------------------------------------------------------------===//

#ifndef KLEE_SEARCHER_H
#define KLEE_SEARCHER_H

#include "klee/Internal/System/Time.h"

#include "llvm/Support/raw_ostream.h"

#include <map>
#include <queue>
#include <set>
#include <unordered_map>
#include <vector>
#include <functional>

namespace llvm {
  class BasicBlock;
  class Function;
  class Instruction;
  class raw_ostream;
}

namespace klee {
  template<class T> class DiscretePDF;
  class ExecutionState;
  class Executor;

  class Searcher {
  public:
    virtual ~Searcher();

    virtual ExecutionState &selectState() = 0;

    virtual void update(ExecutionState *current,
                        const std::vector<ExecutionState *> &addedStates,
                        const std::vector<ExecutionState *> &removedStates) = 0;
    virtual std::vector<ExecutionState *> selectForDelition(int size) { return {}; };

    virtual bool empty() = 0;
    virtual int getSize() {return 0; }

    // prints name of searcher as a klee_message()
    // TODO: could probably make prettier or more flexible
    virtual void printName(llvm::raw_ostream &os) {
      os << "<unnamed searcher>\n";
    }

    // pgbovine - to be called when a searcher gets activated and
    // deactivated, say, by a higher-level searcher; most searchers
    // don't need this functionality, so don't have to override.
    virtual void activate() {}
    virtual void deactivate() {}

    // utility functions

    void addState(ExecutionState *es, ExecutionState *current = 0) {
      std::vector<ExecutionState *> tmp;
      tmp.push_back(es);
      update(current, tmp, std::vector<ExecutionState *>());
    }

    void removeState(ExecutionState *es, ExecutionState *current = 0) {
      std::vector<ExecutionState *> tmp;
      tmp.push_back(es);
      update(current, std::vector<ExecutionState *>(), tmp);
    }

    enum CoreSearchType {
      DFS,
      BFS,
      RandomState,
      RandomPath,
      ZESTI,
      NURS_CovNew,
      NURS_MD2U,
      NURS_Depth,
      NURS_RP,
      NURS_ICnt,
      NURS_CPICnt,
      NURS_QC
    };
  };

  class DFSSearcher : public Searcher {
    std::vector<ExecutionState*> states;

  public:
    ExecutionState &selectState();
    void update(ExecutionState *current,
                const std::vector<ExecutionState *> &addedStates,
                const std::vector<ExecutionState *> &removedStates);
    bool empty() { return states.empty(); }
    void printName(llvm::raw_ostream &os) {
      os << "DFSSearcher\n";
    }
  };

  class EmptySearcher : public Searcher {
  public:
    ExecutionState &selectState() { assert(0 && "Empty searcher is always empty");}
    void update(ExecutionState *current,
                const std::vector<ExecutionState *> &addedStates,
                const std::vector<ExecutionState *> &removedStates) {}
    bool empty() { return true; }
    void printName(llvm::raw_ostream &os) {
      os << "EmptySearcher\n";
    }
  };

  class SwappingSearcher : public Searcher {
    Searcher* searchers[2];
    unsigned currentSearcher = 0;
    std::function<void(void)> swapCallback;

  public:
    SwappingSearcher(Searcher* s1, Searcher* s2, std::function<void(void)> _cb) : swapCallback(_cb) {
        searchers[0] = s1;
        searchers[1] = s2;
    }
    virtual ~SwappingSearcher() { delete searchers[0]; delete searchers[1]; }
    ExecutionState &selectState() { return searchers[currentSearcher]->selectState(); }
    void update(ExecutionState *current,
                const std::vector<ExecutionState *> &addedStates,
                const std::vector<ExecutionState *> &removedStates) {
      for (int i = currentSearcher; i < 2; i++)
        (searchers[i])->update(current, addedStates, removedStates);
    }
    bool empty() { 
        auto ret = searchers[currentSearcher]->empty();
        if(currentSearcher == 0 && ret) {
            currentSearcher++;
            swapCallback();
            //TODO: add callback to disable pending mode
            return empty();
        }
        return ret;
    }
    void printName(llvm::raw_ostream &os) {
      os << "<SwappingSearcher>\n";
      for (int i = 0; i < 2; i++)
        (searchers[i])->printName(os);
      os << "</SwappingSearcher>\n";

    }
  };

  class BFSSearcher : public Searcher {
    std::deque<ExecutionState*> states;

  public:
    ExecutionState &selectState();
    void update(ExecutionState *current,
                const std::vector<ExecutionState *> &addedStates,
                const std::vector<ExecutionState *> &removedStates);
    bool empty() { return states.empty(); }
    void printName(llvm::raw_ostream &os) {
      os << " BFSSearcher size: " << states.size() << "\n";
    }
  };

  class RandomSearcher : public Searcher {
    std::vector<ExecutionState*> states;

  public:
    ExecutionState &selectState();
    void update(ExecutionState *current,
                const std::vector<ExecutionState *> &addedStates,
                const std::vector<ExecutionState *> &removedStates);
    bool empty() { return states.empty(); }
    void printName(llvm::raw_ostream &os) {
      os << "RandomSearcher\n";
    }
  };

  class WeightedRandomSearcher : public Searcher {
  public:
    enum WeightType {
      Depth,
      RP,
      QueryCost,
      InstCount,
      CPInstCount,
      MinDistToUncovered,
      CoveringNew
    };

  private:
    DiscretePDF<ExecutionState*> *states;
    WeightType type;
    bool updateWeights;
    unsigned size = 0;
    
    double getWeight(ExecutionState*);

  public:
    WeightedRandomSearcher(WeightType type);
    ~WeightedRandomSearcher();

    ExecutionState &selectState();
    std::vector<ExecutionState *> selectForDelition(int size);
    void update(ExecutionState *current,
                const std::vector<ExecutionState *> &addedStates,
                const std::vector<ExecutionState *> &removedStates);
    bool empty();
    void printName(llvm::raw_ostream &os) {
      os << "(NURS:";
      switch(type) {
      case Depth              : os << "D "; break; 
      case RP                 : os << "R "; break; 
      case QueryCost          : os << "Q "; break; 
      case InstCount          : os << "I "; break; 
      case CPInstCount        : os << "C "; break; 
      case MinDistToUncovered : os << "M "; break; 
      case CoveringNew        : os << "C "; break; 
      default                 : os << "<unknown type>\n"; return;
      } 
      os << size << ")";
    }
     int getSize() {return size; }

  };

  class RandomPathSearcher : public Searcher {
    Executor &executor;
    int size = 0;
    static int numRPSearchers;
    const uint8_t idBitMask;

  public:
    RandomPathSearcher(Executor &_executor);
    ~RandomPathSearcher();

    ExecutionState &selectState();
    std::vector<ExecutionState *> selectForDelition(int size);

    void update(ExecutionState *current,
                const std::vector<ExecutionState *> &addedStates,
                const std::vector<ExecutionState *> &removedStates);
    bool empty();
    void printName(llvm::raw_ostream &os) {

      os << "(RP " << size << ")";
    }
     int getSize() {return size; }
  };

  class MergeHandler;
  class MergingSearcher : public Searcher {
    friend class MergeHandler;

    private:

    Executor &executor;
    Searcher *baseSearcher;

    public:
    MergingSearcher(Executor &executor, Searcher *baseSearcher);
    ~MergingSearcher();

    ExecutionState &selectState();

    void update(ExecutionState *current,
                const std::vector<ExecutionState *> &addedStates,
                const std::vector<ExecutionState *> &removedStates) {
      baseSearcher->update(current, addedStates, removedStates);
    }
    bool empty() { return baseSearcher->empty(); }
    void printName(llvm::raw_ostream &os) {
      os << "MergingSearcher\n";
    }
  };
  class BatchingSearcher : public Searcher {
    Searcher *baseSearcher;
    time::Span timeBudget;
    unsigned instructionBudget;

    ExecutionState *lastState;
    time::Point lastStartTime;
    unsigned lastStartInstructions;

  public:
    BatchingSearcher(Searcher *baseSearcher, 
                     time::Span _timeBudget,
                     unsigned _instructionBudget);
    ~BatchingSearcher();

    ExecutionState &selectState();
    void update(ExecutionState *current,
                const std::vector<ExecutionState *> &addedStates,
                const std::vector<ExecutionState *> &removedStates);
    bool empty() { return baseSearcher->empty(); }
    void printName(llvm::raw_ostream &os) {
      os << "<BatchingSearcher> timeBudget: " << timeBudget
         << ", instructionBudget: " << instructionBudget
         << ", baseSearcher:\n";
      baseSearcher->printName(os);
      os << "</BatchingSearcher>\n";
    }
  };

  class PendingSearcher : public Searcher {
    Searcher *baseNormalSearcher;
    Searcher *basePendingSearcher;
    Executor* exec;
    time::Span maxReviveTime;

  public:
    PendingSearcher(Searcher *baseNormalSearcher, Searcher *basePendingSearcher, Executor* exec);
    ~PendingSearcher();

    ExecutionState &selectState();
    std::vector<ExecutionState *> selectForDelition(int size);
    void update(ExecutionState *current,
                const std::vector<ExecutionState *> &addedStates,
                const std::vector<ExecutionState *> &removedStates);
    bool empty();
    void printName(llvm::raw_ostream &os) {
      os << "<PendingSearcher>\n";
      baseNormalSearcher->printName(os);
      basePendingSearcher->printName(os);
      os << "</PendingSearcher>\n";
    }
  };

  class ZESTIPendingSearcher : public Searcher {
    Executor* exec;
    Searcher *normalSearcher;
    //Base depth from which we do bounded exploration
    int currentBaseDepth = -1;
    int bound = 0;

    //marks the point where the exploration mode started and it should not accept pending states any more
    bool hasSelectedState = false;
    time::Span maxReviveTime;
    std::unordered_map<const ExecutionState*, int> smallestSensitiveDistance;

    std::vector<ExecutionState*> normalStates;
    std::vector<ExecutionState*> pendingStates;
    std::vector<ExecutionState*> toDelete;

  public:
    ZESTIPendingSearcher(Executor* _exec);
    void computeDistances();

    ExecutionState &selectState();
    void update(ExecutionState *current,
                const std::vector<ExecutionState *> &addedStates,
                const std::vector<ExecutionState *> &removedStates);
    bool empty();
    void printName(llvm::raw_ostream &os) {
      os << "<ZESTIPendingSearcher>\n";
    }
  };


  class IterativeDeepeningTimeSearcher : public Searcher {
    Searcher *baseSearcher;
    time::Point startTime;
    time::Span time;
    std::set<ExecutionState*> pausedStates;

  public:
    IterativeDeepeningTimeSearcher(Searcher *baseSearcher);
    ~IterativeDeepeningTimeSearcher();

    ExecutionState &selectState();
    void update(ExecutionState *current,
                const std::vector<ExecutionState *> &addedStates,
                const std::vector<ExecutionState *> &removedStates);
    bool empty() { return baseSearcher->empty() && pausedStates.empty(); }
    void printName(llvm::raw_ostream &os) {
      os << "IterativeDeepeningTimeSearcher\n";
    }
  };

  class InterleavedSearcher : public Searcher {
    typedef std::vector<Searcher*> searchers_ty;

    searchers_ty searchers;
    unsigned index;

  public:
    explicit InterleavedSearcher(const searchers_ty &_searchers);
    ~InterleavedSearcher();

    ExecutionState &selectState();
    void update(ExecutionState *current,
                const std::vector<ExecutionState *> &addedStates,
                const std::vector<ExecutionState *> &removedStates);
    bool empty() { return searchers[0]->empty(); }
    void printName(llvm::raw_ostream &os) {
      os << "<InterleavedSearcher> containing "
         << searchers.size() << " searchers:\n";
      for (searchers_ty::iterator it = searchers.begin(), ie = searchers.end();
           it != ie; ++it)
        (*it)->printName(os);
      os << "</InterleavedSearcher>\n";
    }
  };

}

#endif /* KLEE_SEARCHER_H */
