//===-- Z3IntSolver.h
//---------------------------------------------------===//
//
//                     The KLEE Symbolic Virtual Machine
//
// This file is distributed under the University of Illinois Open Source
// License. See LICENSE.TXT for details.
//
//===----------------------------------------------------------------------===//

#ifndef KLEE_Z3INTSOLVER_H
#define KLEE_Z3INTSOLVER_H

#include "klee/Solver.h"

namespace klee {
/// Z3IntSolver - A partial int solver based on Z3
class Z3IntSolver : public Solver {
  Solver* coreSolver;
public:
  /// Z3IntSolver - Construct a new Z3IntSolver.
  Z3IntSolver(Solver* s);

  /// Get the query in SMT-LIBv2 format.
  /// \return A C-style string. The caller is responsible for freeing this.
  virtual char *getConstraintLog(const Query &);

  /// setCoreSolverTimeout - Set constraint solver timeout delay to the given
  /// value; 0
  /// is off.
  virtual void setCoreSolverTimeout(time::Span timeout);
};
}

#endif /* KLEE_Z3SOLVER_H */
