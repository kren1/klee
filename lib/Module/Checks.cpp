//===-- Checks.cpp --------------------------------------------------------===//
//
//                     The KLEE Symbolic Virtual Machine
//
// This file is distributed under the University of Illinois Open Source
// License. See LICENSE.TXT for details.
//
//===----------------------------------------------------------------------===//

#include "Passes.h"

#include "klee/Config/Version.h"

#include "KLEEIRMetaData.h"

#include "llvm/Analysis/LoopInfo.h"
#include "llvm/IR/Constants.h"
#include "llvm/IR/DataLayout.h"
#include "llvm/IR/DerivedTypes.h"
#include "llvm/IR/Function.h"
#include "llvm/IR/IRBuilder.h"
#include "llvm/IR/InstrTypes.h"
#include "llvm/IR/Instruction.h"
#include "llvm/IR/Instructions.h"
#include "llvm/IR/IntrinsicInst.h"
#include "llvm/IR/LLVMContext.h"
#include "llvm/IR/Module.h"
#include "llvm/IR/Type.h"
#include "llvm/Pass.h"
#include "llvm/Transforms/Scalar.h"
#include "llvm/Transforms/Utils/BasicBlockUtils.h"
#include "llvm/Analysis/PostDominators.h"

#include<random>

using namespace llvm;
using namespace klee;

cl::opt<double> FaultPct(
    "fault-pct", cl::init(1.),
    cl::desc("percetegae of basic block with faults"
             "(default=1.0 (always))"));

char DivCheckPass::ID;

bool DivCheckPass::runOnModule(Module &M) {
  std::vector<llvm::BinaryOperator *> divInstruction;

  for (auto &F : M) {
    for (auto &BB : F) {
      for (auto &I : BB) {
        auto binOp = dyn_cast<BinaryOperator>(&I);
        if (!binOp)
          continue;

        // find all [s|u][div|rem] instructions
        auto opcode = binOp->getOpcode();
        if (opcode != Instruction::SDiv && opcode != Instruction::UDiv &&
            opcode != Instruction::SRem && opcode != Instruction::URem)
          continue;

        // Check if the operand is constant and not zero, skip in that case.
        const auto &operand = binOp->getOperand(1);
        if (const auto &coOp = dyn_cast<llvm::Constant>(operand)) {
          if (!coOp->isZeroValue())
            continue;
        }

        // Check if the operand is already checked by "klee_div_zero_check"
        if (KleeIRMetaData::hasAnnotation(I, "klee.check.div", "True"))
          continue;
        divInstruction.push_back(binOp);
      }
    }
  }

  // If nothing to do, return
  if (divInstruction.empty())
    return false;

  LLVMContext &ctx = M.getContext();
  KleeIRMetaData md(ctx);
  auto divZeroCheckFunction =
      M.getOrInsertFunction("klee_div_zero_check", Type::getVoidTy(ctx),
                            Type::getInt64Ty(ctx) KLEE_LLVM_GOIF_TERMINATOR);

  for (auto &divInst : divInstruction) {
    llvm::IRBuilder<> Builder(divInst /* Inserts before divInst*/);
    auto denominator =
        Builder.CreateIntCast(divInst->getOperand(1), Type::getInt64Ty(ctx),
                              false, /* sign doesn't matter */
                              "int_cast_to_i64");
    Builder.CreateCall(divZeroCheckFunction, denominator);
    md.addAnnotation(*divInst, "klee.check.div", "True");
  }

  return true;
}

char OvershiftCheckPass::ID;

bool OvershiftCheckPass::runOnModule(Module &M) {
  std::vector<llvm::BinaryOperator *> shiftInstructions;
  for (auto &F : M) {
    for (auto &BB : F) {
      for (auto &I : BB) {
        auto binOp = dyn_cast<BinaryOperator>(&I);
        if (!binOp)
          continue;

        // find all shift instructions
        auto opcode = binOp->getOpcode();
        if (opcode != Instruction::Shl && opcode != Instruction::LShr &&
            opcode != Instruction::AShr)
          continue;

        // Check if the operand is constant and not zero, skip in that case
        auto operand = binOp->getOperand(1);
        if (auto coOp = dyn_cast<llvm::ConstantInt>(operand)) {
          auto typeWidth =
              binOp->getOperand(0)->getType()->getScalarSizeInBits();
          // If the constant shift is positive and smaller,equal the type width,
          // we can ignore this instruction
          if (!coOp->isNegative() && coOp->getZExtValue() < typeWidth)
            continue;
        }

        if (KleeIRMetaData::hasAnnotation(I, "klee.check.shift", "True"))
          continue;

        shiftInstructions.push_back(binOp);
      }
    }
  }

  if (shiftInstructions.empty())
    return false;

  // Retrieve the checker function
  auto &ctx = M.getContext();
  KleeIRMetaData md(ctx);
  auto overshiftCheckFunction = M.getOrInsertFunction(
      "klee_overshift_check", Type::getVoidTy(ctx), Type::getInt64Ty(ctx),
      Type::getInt64Ty(ctx) KLEE_LLVM_GOIF_TERMINATOR);

  for (auto &shiftInst : shiftInstructions) {
    llvm::IRBuilder<> Builder(shiftInst);

    std::vector<llvm::Value *> args;

    // Determine bit width of first operand
    uint64_t bitWidth = shiftInst->getOperand(0)->getType()->getScalarSizeInBits();
    auto bitWidthC = ConstantInt::get(Type::getInt64Ty(ctx), bitWidth, false);
    args.push_back(bitWidthC);

    auto shiftValue =
        Builder.CreateIntCast(shiftInst->getOperand(1), Type::getInt64Ty(ctx),
                              false, /* sign doesn't matter */
                              "int_cast_to_i64");
    args.push_back(shiftValue);

    Builder.CreateCall(overshiftCheckFunction, args);
    md.addAnnotation(*shiftInst, "klee.check.shift", "True");
  }

  return true;
}

char DivFaultPass::ID;

void DivFaultPass::getAnalysisUsage(AnalysisUsage &AU) const {
     AU.setPreservesCFG();
     AU.addRequired<LoopInfoWrapperPass>();
}
bool DivFaultPass::runOnModule(Module &M) {
  LLVMContext &ctx = M.getContext();
  KleeIRMetaData md(ctx);
  APInt location(32, 0);
  std::mt19937 generator(5564354);
  std::bernoulli_distribution distribution(FaultPct);

  auto divFaultFunction = cast<Function>(
      M.getOrInsertFunction("klee_div_fault", Type::getVoidTy(ctx),
                            Type::getInt32Ty(ctx) KLEE_LLVM_GOIF_TERMINATOR));


  for (auto &F : M) {
    if (F.hasName() && F.getName().startswith("klee_"))
      continue;
    if(F.isDeclaration()) continue;
    const DebugLoc* b = nullptr;
    for(auto &BB : F) {
        for(auto &I : BB) {
            if(I.getDebugLoc()) {
                b = &I.getDebugLoc();
                break;
            }
        }
    }

//    LoopInfo &LI = getAnalysis<LoopInfoWrapperPass>(F).getLoopInfo();
    for (auto &BB : F) {
//      if(LI.getLoopFor(&BB) != nullptr) continue; //skip faults in loops
      if (KleeIRMetaData::hasAnnotation(BB.back(), "klee.fault.div", "True"))
        continue;
      if(!distribution(generator)) {
        md.addAnnotation(BB.back(), "klee.fault.div", "True");
        continue;
      }
      location++;
      llvm::IRBuilder<> Builder(&BB.back());
      auto call = Builder.CreateCall(divFaultFunction, Constant::getIntegerValue(Type::getInt32Ty(ctx), location));
      if(b)
        call->setDebugLoc(DebugLoc(b->get()));
      md.addAnnotation(BB.back(), "klee.fault.div", "True");

    }
  }

  return true;
}

char ErrorBBAnnotator::ID;

void ErrorBBAnnotator::getAnalysisUsage(AnalysisUsage &AU) const {
     AU.setPreservesCFG();
     AU.addRequired<PostDominatorTreeWrapperPass>();
}
bool ErrorBBAnnotator::runOnFunction(llvm::Function &F) {
  LLVMContext &ctx = F.getContext();
  KleeIRMetaData md(ctx);
	auto &PDT = getAnalysis<PostDominatorTreeWrapperPass>().getPostDomTree();
  llvm::SmallPtrSet<BasicBlock*, 32> errorBBs;
	for(auto &BB : F) {
  for(auto &I : BB) {
      if(CallInst* ci = dyn_cast<CallInst>(&I)) {
          if(ci->getCalledFunction() && ci->getCalledFunction()->hasName() 
             && ci->getCalledFunction()->getName().startswith("klee_report_error")) {
//                 return true;
                errorBBs.insert(&BB);
             } else if(false && ci->getCalledFunction() && ci->getCalledFunction()->hasName() 
             && ci->getCalledFunction()->getName().startswith("sqlite3ErrorMsg")) {
//                 md.addAnnotation(BB.front(), "klee.sqlerror.block", "True");
                errorBBs.insert(&BB);
                  
                 
             }
    
      }

  }
	}
  
  for(auto BB : errorBBs) {
      llvm::SmallVector<BasicBlock*, 32>  dominated;
      PDT.getDescendants(BB, dominated);
      if(dominated.size() > 1) errs() << "Multiple dominations!\n";
      for(auto domBB : dominated) {
          md.addAnnotation(domBB->front(), "klee.error.block", "True");
      }

  }
  return true;
}
