//===-- Z3IntSolver.cpp -------------------------------------------*- C++ -*-====//
//
//                     The KLEE Symbolic Virtual Machine
//
// This file is distributed under the University of Illinois Open Source
// License. See LICENSE.TXT for details.
//
//===----------------------------------------------------------------------===//
#include "klee/Config/config.h"
#include "klee/Internal/Support/ErrorHandling.h"
#include "klee/Internal/Support/FileHandling.h"
#ifdef ENABLE_Z3
#include "Z3IntSolver.h"
#include "Z3IntBuilder.h"
#include "klee/Constraints.h"
#include "klee/Solver.h"
#include "klee/SolverImpl.h"
#include "klee/util/Assignment.h"
#include "klee/util/ExprUtil.h"
#include "llvm/Support/CommandLine.h"
#include "llvm/Support/raw_ostream.h"

namespace {
// NOTE: Very useful for debugging Z3 behaviour. These files can be given to
// the z3 binary to replay all Z3 API calls using its `-log` option.
llvm::cl::opt<std::string> Z3LogInteractionFile(
    "debug-z3int-log-api-interaction", llvm::cl::init(""),
    llvm::cl::desc("Log API interaction with Z3 to the specified path"));

llvm::cl::opt<std::string> Z3QueryDumpFile(
    "debug-z3int-dump-queries", llvm::cl::init(""),
    llvm::cl::desc("Dump Z3's representation of the query to the specified path"));

llvm::cl::opt<bool> Z3ValidateModels(
    "debug-z3int-validate-models", llvm::cl::init(false),
    llvm::cl::desc("When generating Z3 models validate these against the query"));


llvm::cl::opt<unsigned>
    Z3VerbosityLevel("debug-z3int-verbosity", llvm::cl::init(0),
                     llvm::cl::desc("Z3 verbosity level (default=0)"));
}
bool wasIntQuery = false;

#include "llvm/Support/ErrorHandling.h"

namespace klee {
class IsIntExpr : public ExprVisitor {
  int isIntExpr = true;
  Action no() {
      isIntExpr = false;
      return Action::skipChildren();
    }
  protected:
      
    Action visitXor(const XorExpr& e) override { 
       if(e.right->getWidth() != Expr::Bool) {
        ++stats::xorFail; 
        return no(); 
      }
      return Action::doChildren();
    }
    Action visitShl(const ShlExpr&) override { ++stats::shiftFail; return no(); }
    Action visitLShr(const LShrExpr&) override {++stats::shiftFail;  return no(); }
    Action visitAShr(const AShrExpr&) override {++stats::shiftFail;  return no(); }
    Action visitExtract(const ExtractExpr&) override {++stats::extractFail; return no(); }
    Action visitRead(const ReadExpr& re) override {  
//        if( re.getWidth() != re.updates.root->valueType) return no();

        //ExprVisitor doesn't visit updatelists
        const UpdateNode* un = re.updates.head;
        while(un != nullptr) {
            visit(un->index);
            if(!isIntExpr) return no();
            visit(un->value);
            if(!isIntExpr) return no();
            un = un->next;
        }

        return Action::doChildren();
    }
    Action visitConcat (const ConcatExpr& e) override {
      int readLSB_width;
      const ReadExpr* re = Z3IntBuilder::hasOrderedReads(&e, -1, readLSB_width);
      if(re == nullptr ||
        readLSB_width != re->updates.root->valueType) {
          ++stats::readLSBMissMatchFail;
          return no();
      }
      return visitRead(*re);
    }
 
    Action visitAnd(const AndExpr& e) override { 
      if(e.right->getWidth() != Expr::Bool) {
        ++stats::bitwiseLogicFail;
        return no(); 
      }
      return Action::doChildren();
    }
    Action visitOr(const OrExpr& e) override { 
      if(e.right->getWidth() != Expr::Bool) {
        ++stats::bitwiseLogicFail;
        return no(); 
      }
      return Action::doChildren();
    }
    Action visitNot(const NotExpr& e) override { 
      if(e.expr->getWidth() != Expr::Bool) {
        ++stats::bitwiseLogicFail;
        return no(); 
      }
      return Action::doChildren();
    }     

    explicit operator bool() const { return isIntExpr; }
    IsIntExpr() {}
    void reset() {isIntExpr = true;}
  public:
    static bool isIntQuery(const Query& q) {
        IsIntExpr isInt;
        isInt.visit(q.expr);
        if(!isInt) return false;
        isInt.reset();
        for(auto &e : q.constraints) {
            isInt.visit(e);
            if(!isInt) return false;
            isInt.reset();
        }
        return true;
        
    }
};

class Z3IntSolverImpl : public SolverImpl {
private:
  Z3IntBuilder *builder;
  Solver* coreSolver;
  time::Span timeout;
  SolverRunStatus runStatusCode;
  std::unique_ptr<llvm::raw_fd_ostream> dumpedQueriesFile;
  ::Z3_params solverParameters;
  // Parameter symbols
  ::Z3_symbol timeoutParamStrSymbol;

  bool internalRunSolver(const Query &,
                         const std::vector<const Array *> *objects,
                         std::vector<std::vector<unsigned char> > *values,
                         bool &hasSolution);
  bool validateZ3Model(::Z3_solver &theSolver, ::Z3_model &theModel);

public:
  Z3IntSolverImpl(Solver *s);
  ~Z3IntSolverImpl();

  char *getConstraintLog(const Query &);
  void setCoreSolverTimeout(time::Span _timeout) {
    timeout = _timeout;

    auto timeoutInMilliSeconds = static_cast<unsigned>((timeout.toMicroseconds() / 1000));
    if (!timeoutInMilliSeconds)
      timeoutInMilliSeconds = UINT_MAX;
    Z3_params_set_uint(builder->ctx, solverParameters, timeoutParamStrSymbol,
                       timeoutInMilliSeconds);
  }

  bool computeTruth(const Query &, bool &isValid);
  bool computeValue(const Query &, ref<Expr> &result);
  bool computeInitialValues(const Query &,
                            const std::vector<const Array *> &objects,
                            std::vector<std::vector<unsigned char> > &values,
                            bool &hasSolution);
  SolverRunStatus
  handleSolverResponse(::Z3_solver theSolver, ::Z3_lbool satisfiable,
                       const std::vector<const Array *> *objects,
                       std::vector<std::vector<unsigned char> > *values,
                       bool &hasSolution);
  SolverRunStatus getOperationStatusCode();
};

Z3IntSolverImpl::Z3IntSolverImpl(Solver *s)
    : builder(new Z3IntBuilder(
          /*autoClearConstructCache=*/false,
          /*z3LogInteractionFileArg=*/Z3LogInteractionFile.size() > 0
              ? Z3LogInteractionFile.c_str()
              : NULL)), coreSolver(s),
      runStatusCode(SOLVER_RUN_STATUS_FAILURE) {
  assert(builder && "unable to create Z3IntBuilder");
  solverParameters = Z3_mk_params(builder->ctx);
  Z3_params_inc_ref(builder->ctx, solverParameters);
  timeoutParamStrSymbol = Z3_mk_string_symbol(builder->ctx, "timeout");
  setCoreSolverTimeout(timeout);

  if (!Z3QueryDumpFile.empty()) {
    std::string error;
    dumpedQueriesFile = klee_open_output_file(Z3QueryDumpFile, error);
    if (!dumpedQueriesFile) {
      klee_error("Error creating file for dumping Z3 queries: %s",
                 error.c_str());
    }
    klee_message("Dumping Z3 queries to \"%s\"", Z3QueryDumpFile.c_str());
  }

  // Set verbosity
  if (Z3VerbosityLevel > 0) {
    std::string underlyingString;
    llvm::raw_string_ostream ss(underlyingString);
    ss << Z3VerbosityLevel;
    ss.flush();
    Z3_global_param_set("verbose", underlyingString.c_str());
  }
}

Z3IntSolverImpl::~Z3IntSolverImpl() {
  Z3_params_dec_ref(builder->ctx, solverParameters);
  delete builder;
}

Z3IntSolver::Z3IntSolver(Solver *s) : Solver(new Z3IntSolverImpl(s)), coreSolver(s) {}

char *Z3IntSolver::getConstraintLog(const Query &query) {
  return impl->getConstraintLog(query);
}

void Z3IntSolver::setCoreSolverTimeout(time::Span timeout) {
  impl->setCoreSolverTimeout(timeout);
}

char *Z3IntSolverImpl::getConstraintLog(const Query &query) {
  std::vector<Z3ASTHandle> assumptions;
  // We use a different builder here because we don't want to interfere
  // with the solver's builder because it may change the solver builder's
  // cache.
  // NOTE: The builder does not set `z3LogInteractionFile` to avoid conflicting
  // with whatever the solver's builder is set to do.
  Z3IntBuilder temp_builder(/*autoClearConstructCache=*/false,
                         /*z3LogInteractionFile=*/NULL);
  ConstantArrayFinder constant_arrays_in_query;
  for (auto const &constraint : query.constraints) {
    assumptions.push_back(temp_builder.construct(constraint));
    constant_arrays_in_query.visit(constraint);
  }

  // KLEE Queries are validity queries i.e.
  // ∀ X Constraints(X) → query(X)
  // but Z3 works in terms of satisfiability so instead we ask the
  // the negation of the equivalent i.e.
  // ∃ X Constraints(X) ∧ ¬ query(X)
  Z3ASTHandle formula = Z3ASTHandle(
      Z3_mk_not(temp_builder.ctx, temp_builder.construct(query.expr)),
      temp_builder.ctx);
  constant_arrays_in_query.visit(query.expr);

  for (auto const &constant_array : constant_arrays_in_query.results) {
    assert(builder->constant_array_assertions.count(constant_array) == 1 &&
           "Constant array found in query, but not handled by Z3IntBuilder");
    for (auto const &arrayIndexValueExpr :
         builder->constant_array_assertions[constant_array]) {
      assumptions.push_back(arrayIndexValueExpr);
    }
  }
  for (auto const &sym_array : constant_arrays_in_query.symbolicArrays) {
    assert(builder->constant_array_assertions.count(sym_array) == 1 &&
           "symbolic array found in query, but not handled by Z3IntBuilder");
    for (auto const &arrayIndexValueExpr :
         builder->constant_array_assertions[sym_array]) {
      assumptions.push_back(arrayIndexValueExpr);
    }
  }
  ::Z3_ast *assumptionsArray = NULL;
  int numAssumptions = assumptions.size();
  if (numAssumptions) {
    assumptionsArray = (::Z3_ast *)malloc(sizeof(::Z3_ast) * numAssumptions);
    for (int index = 0; index < numAssumptions; ++index) {
      assumptionsArray[index] = (::Z3_ast)assumptions[index];
    }
  }

  ::Z3_string result = Z3_benchmark_to_smtlib_string(
      temp_builder.ctx,
      /*name=*/"Emited by klee::Z3IntSolverImpl::getConstraintLog()",
      /*logic=*/"",
      /*status=*/"unknown",
      /*attributes=*/"",
      /*num_assumptions=*/numAssumptions,
      /*assumptions=*/assumptionsArray,
      /*formula=*/formula);

  if (numAssumptions)
    free(assumptionsArray);

  // We need to trigger a dereference before the `temp_builder` gets destroyed.
  // We do this indirectly by emptying `assumptions` and assigning to
  // `formula`.
  assumptions.clear();
  formula = Z3ASTHandle(NULL, temp_builder.ctx);
  // Client is responsible for freeing the returned C-string
  return strdup(result);
}

bool Z3IntSolverImpl::computeTruth(const Query &query, bool &isValid) {
  bool hasSolution = false; // to remove compiler warning
  bool status =
      internalRunSolver(query, /*objects=*/NULL, /*values=*/NULL, hasSolution);
  isValid = !hasSolution;
  return status;
}

bool Z3IntSolverImpl::computeValue(const Query &query, ref<Expr> &result) {
  std::vector<const Array *> objects;
  std::vector<std::vector<unsigned char> > values;
  bool hasSolution;

  // Find the object used in the expression, and compute an assignment
  // for them.
  findSymbolicObjects(query.expr, objects);
  if (!computeInitialValues(query.withFalse(), objects, values, hasSolution))
    return false;
  assert(hasSolution && "state has invalid constraint set");

  // Evaluate the expression with the computed assignment.
  Assignment a(objects, values);
  result = a.evaluate(query.expr);

  return true;
}

bool Z3IntSolverImpl::computeInitialValues(
    const Query &query, const std::vector<const Array *> &objects,
    std::vector<std::vector<unsigned char> > &values, bool &hasSolution) {
  return internalRunSolver(query, &objects, &values, hasSolution);
}

bool Z3IntSolverImpl::internalRunSolver(
    const Query &query, const std::vector<const Array *> *objects,
    std::vector<std::vector<unsigned char> > *values, bool &hasSolution) {

  if(!IsIntExpr::isIntQuery(query)) {
//    llvm::errs() << "Int failed, falling back to coreSolver\n";
    wasIntQuery = false;
    ++stats::nonintQueries;
    return coreSolver->impl->computeInitialValues(query, *objects, *values, hasSolution);
  } else {
    wasIntQuery = true;
    //llvm::errs() << "Using int solver! :)\n";
    ++stats::intQueries;
  }

  TimerStatIncrementer t(stats::queryTime);
  // NOTE: Z3 will switch to using a slower solver internally if push/pop are
  // used so for now it is likely that creating a new solver each time is the
  // right way to go until Z3 changes its behaviour.
  //
  // TODO: Investigate using a custom tactic as described in
  // https://github.com/klee/klee/issues/653
  Z3_solver theSolver = Z3_mk_solver(builder->ctx);
  Z3_solver_inc_ref(builder->ctx, theSolver);
  Z3_solver_set_params(builder->ctx, theSolver, solverParameters);

  runStatusCode = SOLVER_RUN_STATUS_FAILURE;

  ConstantArrayFinder constant_arrays_in_query;
  for (auto const &constraint : query.constraints) {
    Z3_solver_assert(builder->ctx, theSolver, builder->construct(constraint));
    constant_arrays_in_query.visit(constraint);
  }
  ++stats::queries;
  if (objects)
    ++stats::queryCounterexamples;

  Z3ASTHandle z3QueryExpr =
      Z3ASTHandle(builder->construct(query.expr), builder->ctx);
  constant_arrays_in_query.visit(query.expr);

  for (auto const &constant_array : constant_arrays_in_query.results) {
    assert(builder->constant_array_assertions.count(constant_array) == 1 &&
           "Constant array found in query, but not handled by Z3IntBuilder");
    for (auto const &arrayIndexValueExpr :
         builder->constant_array_assertions[constant_array]) {
      Z3_solver_assert(builder->ctx, theSolver, arrayIndexValueExpr);
    }
  }

  for (auto const &sym_array : constant_arrays_in_query.symbolicArrays) {
    assert(builder->constant_array_assertions.count(sym_array) == 1 &&
           "symbolic array found in query, but not handled by Z3IntBuilder");
    for (auto const &arrayIndexValueExpr :
         builder->constant_array_assertions[sym_array]) {
      Z3_solver_assert(builder->ctx, theSolver, arrayIndexValueExpr);
    }
  }
  // KLEE Queries are validity queries i.e.
  // ∀ X Constraints(X) → query(X)
  // but Z3 works in terms of satisfiability so instead we ask the
  // negation of the equivalent i.e.
  // ∃ X Constraints(X) ∧ ¬ query(X)
  Z3_solver_assert(
      builder->ctx, theSolver,
      Z3ASTHandle(Z3_mk_not(builder->ctx, z3QueryExpr), builder->ctx));

  if (dumpedQueriesFile) {
    *dumpedQueriesFile << "; start Z3 query\n";
    *dumpedQueriesFile << Z3_solver_to_string(builder->ctx, theSolver);
    *dumpedQueriesFile << "(check-sat)\n";
    *dumpedQueriesFile << "(reset)\n";
    *dumpedQueriesFile << "; end Z3 query\n\n";
    dumpedQueriesFile->flush();
  }

  ::Z3_lbool satisfiable = Z3_solver_check(builder->ctx, theSolver);
  runStatusCode = handleSolverResponse(theSolver, satisfiable, objects, values,
                                       hasSolution);
  if(runStatusCode == SolverImpl::SOLVER_RUN_STATUS_FAILURE) {
     return coreSolver->impl->computeInitialValues(query, *objects, *values, hasSolution);
  }

  Z3_solver_dec_ref(builder->ctx, theSolver);
  // Clear the builder's cache to prevent memory usage exploding.
  // By using ``autoClearConstructCache=false`` and clearning now
  // we allow Z3_ast expressions to be shared from an entire
  // ``Query`` rather than only sharing within a single call to
  // ``builder->construct()``.
  builder->clearConstructCache();

  if (runStatusCode == SolverImpl::SOLVER_RUN_STATUS_SUCCESS_SOLVABLE ||
      runStatusCode == SolverImpl::SOLVER_RUN_STATUS_SUCCESS_UNSOLVABLE) {
    if (hasSolution) {
      ++stats::queriesInvalid;
    } else {
      ++stats::queriesValid;
    }
    return true; // success
  }
  return false; // failed
}

SolverImpl::SolverRunStatus Z3IntSolverImpl::handleSolverResponse(
    ::Z3_solver theSolver, ::Z3_lbool satisfiable,
    const std::vector<const Array *> *objects,
    std::vector<std::vector<unsigned char> > *values, bool &hasSolution) {
  switch (satisfiable) {
  case Z3_L_TRUE: {
    hasSolution = true;
    if (!objects) {
      // No assignment is needed
      assert(values == NULL);
      return SolverImpl::SOLVER_RUN_STATUS_SUCCESS_SOLVABLE;
    }
    assert(values && "values cannot be nullptr");
    ::Z3_model theModel = Z3_solver_get_model(builder->ctx, theSolver);
    assert(theModel && "Failed to retrieve model");
    Z3_model_inc_ref(builder->ctx, theModel);
    values->reserve(objects->size());
    for (std::vector<const Array *>::const_iterator it = objects->begin(),
                                                    ie = objects->end();
         it != ie; ++it) {
      const Array *array = *it;
//      llvm::errs() << "Array type: " << array->valueType << "\n";
      assert(array->valueType != Expr::InvalidWidth && "Int solver needs value type");
      unsigned byteStride = array->valueType / 8;
      std::vector<unsigned char> data;

      data.reserve(array->size);
      for (unsigned offset = 0; offset < (array->size / byteStride); offset++) {
        // We can't use Z3ASTHandle here so have to do ref counting manually
        ::Z3_ast arrayElementExpr;
        Z3ASTHandle initial_read = builder->getInitialRead(array, offset);

        __attribute__((unused))
        bool successfulEval =
            Z3_model_eval(builder->ctx, theModel, initial_read,
                          /*model_completion=*/Z3_TRUE, &arrayElementExpr);
        assert(successfulEval && "Failed to evaluate model");
        Z3_inc_ref(builder->ctx, arrayElementExpr);
        assert(Z3_get_ast_kind(builder->ctx, arrayElementExpr) ==
                   Z3_NUMERAL_AST &&
               "Evaluated expression has wrong sort");

        int64_t arrayElementValue = 0;
        __attribute__((unused))
        bool successGet = Z3_get_numeral_int64(builder->ctx, arrayElementExpr,
                                             &arrayElementValue);
        
        assert(successGet && "failed to get value back");
        switch(byteStride) {
            case 1: 
                  assert(arrayElementValue <= std::numeric_limits<std::int8_t>::max());
                  assert(arrayElementValue >= std::numeric_limits<std::int8_t>::min());
                  break;
            case 2:
                  assert(arrayElementValue <= std::numeric_limits<std::int16_t>::max()); 
                  assert(arrayElementValue >= std::numeric_limits<std::int16_t>::min()); 
                  break;
            case 4: 
                  assert(arrayElementValue <= std::numeric_limits<std::int32_t>::max());
                  assert(arrayElementValue >= std::numeric_limits<std::int32_t>::min());
                  break;
            case 8:
                  assert(arrayElementValue <= std::numeric_limits<std::int64_t>::max());
                  assert(arrayElementValue >= std::numeric_limits<std::int64_t>::min());
                  break;
            default: assert(0 && "Unknown byte stride");
        }
        uint8_t *p = (uint8_t*)&arrayElementValue;
        //llvm::errs() << array->name << "->" << arrayElementValue << " bs:" << byteStride  << " p:" << (int)*p <<"\n";
        for(unsigned j = 0; j < byteStride; j++)
          data.push_back(p[j]);
        Z3_dec_ref(builder->ctx, arrayElementExpr);
      }
      values->push_back(data);
    }

    // Validate the model if requested
    if (Z3ValidateModels) {
      bool success = validateZ3Model(theSolver, theModel);
      if (!success)
        abort();
    }

    Z3_model_dec_ref(builder->ctx, theModel);
    return SolverImpl::SOLVER_RUN_STATUS_SUCCESS_SOLVABLE;
  }
  case Z3_L_FALSE:
    hasSolution = false;
    return SolverImpl::SOLVER_RUN_STATUS_SUCCESS_UNSOLVABLE;
  case Z3_L_UNDEF: {
    ::Z3_string reason =
        ::Z3_solver_get_reason_unknown(builder->ctx, theSolver);
    if (strcmp(reason, "timeout") == 0 || strcmp(reason, "canceled") == 0 ||
        strcmp(reason, "(resource limits reached)") == 0) {
      return SolverImpl::SOLVER_RUN_STATUS_TIMEOUT;
    }
    if (strcmp(reason, "unknown") == 0) {
      return SolverImpl::SOLVER_RUN_STATUS_FAILURE;
    }
    klee_warning("Unexpected solver failure. Reason is \"%s,\"\n", reason);
    abort();
  }
  default:
    llvm_unreachable("unhandled Z3 result");
  }
}

bool Z3IntSolverImpl::validateZ3Model(::Z3_solver &theSolver, ::Z3_model &theModel) {
  bool success = true;
  ::Z3_ast_vector constraints =
      Z3_solver_get_assertions(builder->ctx, theSolver);
  Z3_ast_vector_inc_ref(builder->ctx, constraints);

  unsigned size = Z3_ast_vector_size(builder->ctx, constraints);

  for (unsigned index = 0; index < size; ++index) {
    Z3ASTHandle constraint = Z3ASTHandle(
        Z3_ast_vector_get(builder->ctx, constraints, index), builder->ctx);

    ::Z3_ast rawEvaluatedExpr;
    __attribute__((unused))
    bool successfulEval =
        Z3_model_eval(builder->ctx, theModel, constraint,
                      /*model_completion=*/Z3_TRUE, &rawEvaluatedExpr);
    assert(successfulEval && "Failed to evaluate model");

    // Use handle to do ref-counting.
    Z3ASTHandle evaluatedExpr(rawEvaluatedExpr, builder->ctx);

    Z3SortHandle sort =
        Z3SortHandle(Z3_get_sort(builder->ctx, evaluatedExpr), builder->ctx);
    assert(Z3_get_sort_kind(builder->ctx, sort) == Z3_BOOL_SORT &&
           "Evaluated expression has wrong sort");

    Z3_lbool evaluatedValue =
        Z3_get_bool_value(builder->ctx, evaluatedExpr);
    if (evaluatedValue != Z3_L_TRUE) {
      llvm::errs() << "Validating model failed:\n"
                   << "The expression:\n";
      constraint.dump();
      llvm::errs() << "evaluated to \n";
      evaluatedExpr.dump();
      llvm::errs() << "But should be true\n";
      success = false;
    }
  }

  if (!success) {
    llvm::errs() << "Solver state:\n" << Z3_solver_to_string(builder->ctx, theSolver) << "\n";
    llvm::errs() << "Model:\n" << Z3_model_to_string(builder->ctx, theModel) << "\n";
  }

  Z3_ast_vector_dec_ref(builder->ctx, constraints);
  return success;
}

SolverImpl::SolverRunStatus Z3IntSolverImpl::getOperationStatusCode() {
  return runStatusCode;
}
}
#endif // ENABLE_Z3
