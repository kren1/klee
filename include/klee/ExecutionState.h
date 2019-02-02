//===-- ExecutionState.h ----------------------------------------*- C++ -*-===//
//
//                     The KLEE Symbolic Virtual Machine
//
// This file is distributed under the University of Illinois Open Source
// License. See LICENSE.TXT for details.
//
//===----------------------------------------------------------------------===//

#ifndef KLEE_EXECUTIONSTATE_H
#define KLEE_EXECUTIONSTATE_H

#include "klee/Constraints.h"
#include "klee/Expr.h"
#include "klee/Internal/ADT/TreeStream.h"
#include "klee/Internal/System/Time.h"
#include "klee/Internal/Module/Cell.h"
#include "klee/MergeHandler.h"

// FIXME: We do not want to be exposing these? :(
#include "../../lib/Core/AddressSpace.h"
#include "klee/Internal/Module/KInstIterator.h"
#include "klee/Solver.h"

#include <map>
#include <memory>
#include <set>
#include <vector>

namespace klee {
class Array;
class CallPathNode;
struct KFunction;
struct KInstruction;
class MemoryObject;
class PTreeNode;
struct InstructionInfo;

llvm::raw_ostream &operator<<(llvm::raw_ostream &os, const MemoryMap &mm);

struct StackFrame {
  KInstIterator caller;
  KFunction *kf;
  CallPathNode *callPathNode;

  std::vector<const MemoryObject *> allocas;

private:
  std::vector<Cell> locals;
  std::vector<size_t> local_index;
  std::vector<Cell> local_value;

public:
  Cell &getCell(size_t index) {
    if (!locals.empty())
      return locals[index];

    auto pos = std::lower_bound(local_index.begin(), local_index.end(), index);
    auto diff = pos - local_index.begin();
    if (pos != local_index.end() && *pos == index)
      return local_value[diff];

    // insert new value;
    local_index.insert(pos, index);
    local_value.insert(local_value.begin() + diff, Cell());

    return local_value[diff];
  }

  const Cell &getCell(size_t index) const {
    if (!locals.empty())
      return locals[index];

    auto pos = std::lower_bound(local_index.begin(), local_index.end(), index);
    auto diff = pos - local_index.begin();
    if (pos != local_index.end() && *pos == index)
      return local_value[diff];

    assert(0);
    return local_value[diff];
  }

  /// Minimum distance to an uncovered instruction once the function
  /// returns. This is not a good place for this but is used to
  /// quickly compute the context sensitive minimum distance to an
  /// uncovered instruction. This value is updated by the StatsTracker
  /// periodically.
  unsigned minDistToUncoveredOnReturn;

  // For vararg functions: arguments not passed via parameter are
  // stored (packed tightly) in a local (alloca) memory object. This
  // is set up to match the way the front-end generates vaarg code (it
  // does not pass vaarg through as expected). VACopy is lowered inside
  // of intrinsic lowering.
  MemoryObject *varargs;

  StackFrame(KInstIterator caller, KFunction *kf);
};

struct Symbol_t;
struct Symbol_t {
  ref<const MemoryObject> mo;
  const Array *array;

  ref<Symbol_t> next;
  ReferenceCounter __refCount;
  Symbol_t(ref<const MemoryObject> m, const Array *a, ref<Symbol_t> s)
      : mo(std::move(m)), array(a), next(std::move(next)) {}
};

/// @brief ExecutionState representing a path under exploration
class ExecutionState {
public:
  typedef std::vector<StackFrame> stack_ty;

private:
  std::map<std::string, std::string> fnAliases;

public:
  // Execution - Control Flow specific

  /// @brief Pointer to instruction to be executed after the current
  /// instruction
  KInstIterator pc;

  /// @brief Pointer to instruction which is currently executed
  KInstIterator prevPC;

  /// @brief Stack representing the current instruction stream
  stack_ty stack;

  /// @brief Remember from which Basic Block control flow arrived
  /// (i.e. to select the right phi values)
  unsigned incomingBBIndex;

  // Overall state of the state - Data specific

  /// @brief Address space used by this state (e.g. Global and Heap)
  AddressSpace addressSpace;

  /// @brief Constraints collected so far
  SimpleConstraintSet constraints;

  /// Statistics and information

  /// Meta-data for the solver for state-specific queries
  SolverQueryMetaData solverMetaData;

  /// @brief Weight assigned for importance of this state.  Can be
  /// used for searchers to decide what paths to explore
  double weight;

  /// @brief Exploration depth, i.e., number of times KLEE branched for this state
  unsigned depth;

  /// @brief History of complete path: represents branches taken to
  /// reach/create this state (both concrete and symbolic)
  TreeOStream pathOS;

  /// @brief History of symbolic path: represents symbolic branches
  /// taken to reach/create this state
  TreeOStream symPathOS;

  /// @brief Counts how many instructions were executed since the last new
  /// instruction was covered.
  unsigned instsSinceCovNew;

  /// @brief Whether a new instruction was covered in this state
  bool coveredNew;

  /// @brief Disables forking for this state. Set by user code
  bool forkDisabled;

  /// @brief Set containing which lines in which files are covered by this state
  std::map<const std::string *, std::set<unsigned> > coveredLines;

  /// @brief Pointer to the process tree of the current state
  PTreeNode *ptreeNode;

  /// @brief Ordered list of symbolics: used to generate test cases.
  ref<Symbol_t> symbolics;

  std::string getFnAlias(std::string fn);
  void addFnAlias(std::string old_fn, std::string new_fn);
  void removeFnAlias(std::string fn);

  // The objects handling the klee_open_merge calls this state ran through
  std::vector<ref<MergeHandler> > openMergeStack;

  // The numbers of times this state has run through Executor::stepInstruction
  std::uint64_t steppedInstructions;

  uint64_t arrayCntr = 0;

private:
  ExecutionState(const ExecutionState &state);

public:
  ExecutionState() = delete;
  ExecutionState &operator=(const ExecutionState &) = delete;

  ExecutionState(KFunction *kf);

  ~ExecutionState();

  ExecutionState *branch();

  void pushFrame(KInstIterator caller, KFunction *kf);
  void popFrame();

  void addSymbolic(const MemoryObject *mo, const Array *array);
  void addConstraint(ref<Expr> e);

  bool merge(const ExecutionState &b);
  void dumpStack(llvm::raw_ostream &out) const;
};
}

#endif
