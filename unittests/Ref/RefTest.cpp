//===-- RefTest.cpp ---------------------------------------------*- C++ -*-===//
//
//                     The KLEE Symbolic Virtual Machine
//
// This file is distributed under the University of Illinois Open Source
// License. See LICENSE.TXT for details.
//
//===----------------------------------------------------------------------===//

/* Regression test for a bug caused by assigning a ref to itself.
   More details at http://keeda.stanford.edu/pipermail/klee-commits/2012-February/000904.html */

#include "klee/util/Ref.h"
#include "gtest/gtest.h"
#include <iostream>
using klee::ref;

int finished = 0;

struct Expr
{
  /// @brief Required by klee::ref-managed objects
  struct klee::ReferenceCounter __refCount;
  Expr() {
    //std::cout << "Expr(" << this << ") created\n"; 
  }
  ~Expr() { 
    //std::cout << "Expr(" << this << ") destroyed\n"; 
    EXPECT_EQ(finished, 1);
  }
};

TEST(RefTest, SelfAssign) 
{
  struct Expr *r_e = new Expr();
  ref<Expr> r(r_e);
  EXPECT_EQ(r_e->__refCount.refCount, 1u);
  r = r;
  EXPECT_EQ(r_e->__refCount.refCount, 1u);
  finished = 1;
}
