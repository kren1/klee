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

#include "ExecutionState.h"
#include "PTree.h"
#include "klee/ADT/RNG.h"
#include "klee/Support/Debug.h"
#include "klee/System/Time.h"

#include "llvm/Support/CommandLine.h"
#include "llvm/Support/raw_ostream.h"

#include <map>
#include <queue>
#include <set>
#include <unordered_map>
#include <vector>

namespace llvm {
  class BasicBlock;
  class Function;
  class Instruction;
  class raw_ostream;
}

namespace klee {
  template<class T, class Comparator> class DiscretePDF;
  class ExecutionState;
  class Executor;

  /// A Searcher implements an exploration strategy for the Executor by selecting
  /// states for further exploration using different strategies or heuristics.
  class Searcher {
  public:
    virtual ~Searcher() = default;

    /// Selects a state for further exploration.
    /// \return The selected state.
    virtual ExecutionState &selectState() = 0;

    /// Notifies searcher about new or deleted states.
    /// \param current The currently selected state for exploration.
    /// \param addedStates The newly branched states with `current` as common ancestor.
    /// \param removedStates The states that will be terminated.
    virtual void update(ExecutionState *current,
                        const std::vector<ExecutionState *> &addedStates,
                        const std::vector<ExecutionState *> &removedStates) = 0;

    /// \return True if no state left for exploration, False otherwise
    virtual bool empty() = 0;

    /// Prints name of searcher as a `klee_message()`.
    // TODO: could probably made prettier or more flexible
    virtual void printName(llvm::raw_ostream &os) = 0;

    enum CoreSearchType : std::uint8_t {
      DFS,
      BFS,
      RandomState,
      RandomPath,
      NURS_CovNew,
      NURS_MD2U,
      NURS_Depth,
      NURS_RP,
      NURS_ICnt,
      NURS_CPICnt,
      NURS_QC,
      ZESTI
    };
  };

  /// DFSSearcher implements depth-first exploration. All states are kept in
  /// insertion order. The last state is selected for further exploration.
  class DFSSearcher final : public Searcher {
    std::vector<ExecutionState*> states;

  public:
    ExecutionState &selectState() override;
    void update(ExecutionState *current,
                const std::vector<ExecutionState *> &addedStates,
                const std::vector<ExecutionState *> &removedStates) override;
    bool empty() override;
    void printName(llvm::raw_ostream &os) override;
  };

  /// BFSSearcher implements breadth-first exploration. When KLEE branches multiple
  /// times for a single instruction, all new states have the same depth. Keep in
  /// mind that the process tree (PTree) is a binary tree and hence the depth of
  /// a state in that tree and its branch depth during BFS are different.
  class BFSSearcher final : public Searcher {
    std::deque<ExecutionState*> states;

  public:
    ExecutionState &selectState() override;
    void update(ExecutionState *current,
                const std::vector<ExecutionState *> &addedStates,
                const std::vector<ExecutionState *> &removedStates) override;
    bool empty() override;
    void printName(llvm::raw_ostream &os) override;
  };

  /// RandomSearcher picks a state randomly.
  class RandomSearcher final : public Searcher {
    std::vector<ExecutionState*> states;
    RNG &theRNG;

  public:
    explicit RandomSearcher(RNG &rng);
    ExecutionState &selectState() override;
    void update(ExecutionState *current,
                const std::vector<ExecutionState *> &addedStates,
                const std::vector<ExecutionState *> &removedStates) override;
    bool empty() override;
    void printName(llvm::raw_ostream &os) override;
  };

  /// The base class for all weighted searchers. Uses DiscretePDF as underlying
  /// data structure.
  class WeightedRandomSearcher final : public Searcher {
  public:
    enum WeightType : std::uint8_t {
      Depth,
      RP,
      QueryCost,
      InstCount,
      CPInstCount,
      MinDistToUncovered,
      CoveringNew
    };

  private:
    std::unique_ptr<DiscretePDF<ExecutionState*, ExecutionStateIDCompare>> states;
    RNG &theRNG;
    WeightType type;
    bool updateWeights;
    
    double getWeight(ExecutionState*);

  public:
    /// \param type The WeightType that determines the underlying heuristic.
    /// \param RNG A random number generator.
    WeightedRandomSearcher(WeightType type, RNG &rng);
    ~WeightedRandomSearcher() override = default;

    ExecutionState &selectState() override;
    void update(ExecutionState *current,
                const std::vector<ExecutionState *> &addedStates,
                const std::vector<ExecutionState *> &removedStates) override;
    bool empty() override;
    void printName(llvm::raw_ostream &os) override;
  };

  /// RandomPathSearcher performs a random walk of the PTree to select a state.
  /// PTree is a global data structure, however, a searcher can sometimes only
  /// select from a subset of all states (depending on the update calls).
  ///
  /// To support this, RandomPathSearcher has a subgraph view of PTree, in that it
  /// only walks the PTreeNodes that it "owns". Ownership is stored in the
  /// getInt method of the PTreeNodePtr class (which hides it in the pointer itself).
  ///
  /// The current implementation of PTreeNodePtr supports only 3 instances of the
  /// RandomPathSearcher. This is because the current PTreeNodePtr implementation
  /// conforms to C++ and only steals the last 3 alignment bits. This restriction
  /// could be relaxed slightly by an architecture-specific implementation of
  /// PTreeNodePtr that also steals the top bits of the pointer.
  ///
  /// The ownership bits are maintained in the update method.
  class RandomPathSearcher final : public Searcher {
    PTree &processTree;
    RNG &theRNG;

    // Unique bitmask of this searcher
    const uint8_t idBitMask;

  public:
    /// \param processTree The process tree.
    /// \param RNG A random number generator.
    RandomPathSearcher(PTree &processTree, RNG &rng);
    ~RandomPathSearcher() override = default;

    ExecutionState &selectState() override;
    void update(ExecutionState *current,
                const std::vector<ExecutionState *> &addedStates,
                const std::vector<ExecutionState *> &removedStates) override;
    bool empty() override;
    void printName(llvm::raw_ostream &os) override;
  };


  extern llvm::cl::opt<bool> UseIncompleteMerge;
  class MergeHandler;
  class MergingSearcher final : public Searcher {
    friend class MergeHandler;

    private:

    std::unique_ptr<Searcher> baseSearcher;

    /// States that have been paused by the 'pauseState' function
    std::vector<ExecutionState*> pausedStates;

    public:
    /// \param baseSearcher The underlying searcher (takes ownership).
    explicit MergingSearcher(Searcher *baseSearcher);
    ~MergingSearcher() override = default;

    /// ExecutionStates currently paused from scheduling because they are
    /// waiting to be merged in a klee_close_merge instruction
    std::set<ExecutionState *> inCloseMerge;

    /// Keeps track of all currently ongoing merges.
    /// An ongoing merge is a set of states (stored in a MergeHandler object)
    /// which branched from a single state which ran into a klee_open_merge(),
    /// and not all states in the set have reached the corresponding
    /// klee_close_merge() yet.
    std::vector<MergeHandler *> mergeGroups;

    /// Remove state from the searcher chain, while keeping it in the executor.
    /// This is used here to 'freeze' a state while it is waiting for other
    /// states in its merge group to reach the same instruction.
    void pauseState(ExecutionState &state);

    /// Continue a paused state
    void continueState(ExecutionState &state);

    ExecutionState &selectState() override;
    void update(ExecutionState *current,
                const std::vector<ExecutionState *> &addedStates,
                const std::vector<ExecutionState *> &removedStates) override;

    bool empty() override;
    void printName(llvm::raw_ostream &os) override;
  };


  // A simple searcher for pending states. It splits the incoming states into
  // normal and pending states. Pending states are attempted to be revived on
  // calls to update using the fastSolver. If they are succesfully revived they
  // are treated as normal states. The PendingSearcher keeps picking normal
  // states (using baseNormalSearcher) until it can. It then tries to revive
  // pending states (in the order of basePendingSearcher) until it revives one.
  class PendingSearcher : public Searcher {
    // Searcher for normal states
    Searcher *baseNormalSearcher;
    // Searcher for pending states
    Searcher *basePendingSearcher;
    // Used for reviving and terminating infeasible pending states
    Executor &exec;
    // Maximim time spent reviving
    time::Span maxReviveTime;
    // A solver that is used to attempt revives in the update() function. It
    // can be incomplete and should be faster than the normal searcher.
    Solver *fastSolver;

  public:
    PendingSearcher(Searcher *baseNormalSearcher, Searcher *basePendingSearcher,
                    Executor &exec);
    ~PendingSearcher();

    ExecutionState &selectState();
    void update(ExecutionState *current,
                const std::vector<ExecutionState *> &addedStates,
                const std::vector<ExecutionState *> &removedStates);
    bool empty() override;
    void printName(llvm::raw_ostream &os) {
      os << "<PendingSearcher>\n";
      baseNormalSearcher->printName(os);
      basePendingSearcher->printName(os);
      os << "</PendingSearcher>\n";
    }
  };

  /// BatchingSearcher selects a state from an underlying searcher and returns
  /// that state for further exploration for a given time or a given number
  /// of instructions.
  class BatchingSearcher final : public Searcher {
    std::unique_ptr<Searcher> baseSearcher;
    time::Span timeBudget;
    unsigned instructionBudget;

    ExecutionState *lastState {nullptr};
    time::Point lastStartTime;
    unsigned lastStartInstructions;

  public:
    /// \param baseSearcher The underlying searcher (takes ownership).
    /// \param timeBudget Time span a state gets selected before choosing a different one.
    /// \param instructionBudget Number of instructions to re-select a state for.
    BatchingSearcher(Searcher *baseSearcher, time::Span timeBudget, unsigned instructionBudget);
    ~BatchingSearcher() override = default;

    ExecutionState &selectState() override;
    void update(ExecutionState *current,
                const std::vector<ExecutionState *> &addedStates,
                const std::vector<ExecutionState *> &removedStates) override;
    bool empty() override;
    void printName(llvm::raw_ostream &os) override;
  };

  /// IterativeDeepeningTimeSearcher implements time-based deepening. States
  /// are selected from an underlying searcher. When a state reaches its time
  /// limit it is paused (removed from underlying searcher). When the underlying
  /// searcher runs out of states, the time budget is increased and all paused
  /// states are revived (added to underlying searcher).
  class IterativeDeepeningTimeSearcher final : public Searcher {
    std::unique_ptr<Searcher> baseSearcher;
    time::Point startTime;
    time::Span time {time::seconds(1)};
    std::set<ExecutionState*> pausedStates;

  public:
    /// \param baseSearcher The underlying searcher (takes ownership).
    explicit IterativeDeepeningTimeSearcher(Searcher *baseSearcher);
    ~IterativeDeepeningTimeSearcher() override = default;

    ExecutionState &selectState() override;
    void update(ExecutionState *current,
                const std::vector<ExecutionState *> &addedStates,
                const std::vector<ExecutionState *> &removedStates) override;
    bool empty() override;
    void printName(llvm::raw_ostream &os) override;
  };

  /// InterleavedSearcher selects states from a set of searchers in round-robin
  /// manner. It is used for KLEE's default strategy where it switches between
  /// RandomPathSearcher and WeightedRandomSearcher with CoveringNew metric.
  class InterleavedSearcher final : public Searcher {
    std::vector<std::unique_ptr<Searcher>> searchers;
    unsigned index {1};

  public:
    /// \param searchers The underlying searchers (takes ownership).
    explicit InterleavedSearcher(const std::vector<Searcher *> &searchers);
    ~InterleavedSearcher() override = default;

    ExecutionState &selectState() override;
    void update(ExecutionState *current,
                const std::vector<ExecutionState *> &addedStates,
                const std::vector<ExecutionState *> &removedStates) override;
    bool empty() override;
    void printName(llvm::raw_ostream &os) override;
  };

  /// EmptySearcher is always empty
  class EmptySearcher : public Searcher {
  public:
    ExecutionState &selectState() {
      assert(0 && "Empty searcher is always empty");
    }
    void update(ExecutionState *current,
                const std::vector<ExecutionState *> &addedStates,
                const std::vector<ExecutionState *> &removedStates) {}
    bool empty() { return true; }
    void printName(llvm::raw_ostream &os) { os << "EmptySearcher\n"; }
  };

  /// SwappingSearcher has two sub searchers s1 and s2.
  /// It first delegates to s1 untill s1 is empty, then it switches 
  /// Delegation to s2. It calls a callback when the switch occurs.
  /// It delegates the update calls to both searchers at first, but only
  /// to s2 after the switch.
  class SwappingSearcher : public Searcher {
    std::unique_ptr<Searcher> searchers[2];
    unsigned currentSearcher = 0;
    std::function<void(void)> swapCallback;

  public:
    SwappingSearcher(Searcher *s1, Searcher *s2, std::function<void(void)> _cb)
        : swapCallback(_cb) {
      searchers[0] = std::unique_ptr<Searcher>(s1);
      searchers[1] = std::unique_ptr<Searcher>(s2);
    }
    virtual ~SwappingSearcher() override = default;
    ExecutionState &selectState() {
      return searchers[currentSearcher]->selectState();
    }
    void update(ExecutionState *current,
                const std::vector<ExecutionState *> &addedStates,
                const std::vector<ExecutionState *> &removedStates) {
      for (int i = currentSearcher; i < 2; i++)
        (searchers[i])->update(current, addedStates, removedStates);
    }
    bool empty() {
      auto ret = searchers[currentSearcher]->empty();
      if (currentSearcher == 0 && ret) {
        currentSearcher++;
        swapCallback();
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


  ///  ZESTIPendingSearcher implements the ZESTI search "around" sensitive
  ///  instructions using pending states. It assumes a single seed was executed to
  ///  completion with pending states dangling off possible branch points. This
  ///  states can be achived by  
  ///
  ///  SwappingSearcher( 
  ///      PendingSearcher(DFS, Empty),
  ///      ZESTIPendingSearcher())
  ///
  /// Where PendingSearcher(DFS, Empty) just does the seeding and is then
  /// swapped with this searchers. This ZESTI searcher then orders the pending
  /// states it has in the order of decreasing distance from sensitive
  /// instructions, which are then revived in that order. When the revival is
  /// successful, this searcher remebers the depth of that state and starts
  /// killing normal states after that depth + a bound is reached thus achiving
  /// limited bounded exploration.
  class ZESTIPendingSearcher : public Searcher {
    // Used for reviving and terminating states
    Executor &exec; 

    /// Holds the searcher for normal states. As we are doing exhaustive
    /// bounded exploration, this can just be DFS.
    std::unique_ptr<Searcher> normalSearcher;

    /// Depth of the last revived pending states. This is the base depth from
    /// which bounded exploration is performed.
    int currentBaseDepth = -1; 

    /// bound to which we do bounded exploration. States at depth
    /// currentBaseDepth + bound are killed.
    int bound = 0; 

    /// Distance to sensitive instruction multiplier.  The bound above is computed
    /// as ZestiBound * distanceToSensitiveInstruction
    int ZestiBound = 2;

    /// Marks the point where the exploration mode started and it should not
    /// accept pending states any more
    bool hasSelectedState = false;

    /// Maps each (pending) state to the distance to nearest sensitive instruction.
    /// Computed in computeDistances().
    std::unordered_map<const ExecutionState *, int> smallestSensitiveDistance;

    /// Sorted list of all the pending states this searcher knows about.  / The
    /// states are sorted in the order of distance from sensitive instructions by
    /// the computeDistances function.
    std::vector<ExecutionState *> pendingStates;

    /// States that need to be deleted next time it is safe (in selectState)
    std::vector<ExecutionState *> toDelete;

    /// A set of depths where a sensitive instruction occured. Added by
    /// registerSensitiveDepth.
    std::set<std::uint32_t> senstiveDepths;

    /// Computes distances between pending states and their nearest sensitive instructions
    /// and puts in in smallestSensitiveDistance map.
    void computeDistances();

  public:
    ZESTIPendingSearcher(Executor &_exec, int ZestiBound);

    /// Register a depth where sensitive instruction is present.
    void registerSensitiveDepth(std::uint32_t);

    ExecutionState &selectState();
    void update(ExecutionState *current,
                const std::vector<ExecutionState *> &addedStates,
                const std::vector<ExecutionState *> &removedStates);
    bool empty();
    void printName(llvm::raw_ostream &os) { os << "<ZESTIPendingSearcher>\n"; }
  };

} // klee namespace

#endif /* KLEE_SEARCHER_H */
