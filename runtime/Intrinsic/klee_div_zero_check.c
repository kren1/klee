//===-- klee_div_zero_check.c ---------------------------------------------===//
//
//                     The KLEE Symbolic Virtual Machine
//
// This file is distributed under the University of Illinois Open Source
// License. See LICENSE.TXT for details.
//
//===----------------------------------------------------------------------===//

#include "klee/klee.h"

void klee_div_zero_check(long long z) {
  if (z == 0)
    klee_report_error(__FILE__, __LINE__, "divide by zero", "div.err");
}

static int isInit = 0;
static int  symVar;
//void klee_div_fault(int location) {
__attribute__(( always_inline)) void klee_div_fault(int location) {
#pragma clang optimize off
    if(isInit == 0) {
        klee_make_symbolic(&symVar, sizeof(symVar), "divFautlVar");
        isInit = 1;
    }
    if(symVar - location == 0) {
        klee_report_error("fictional", location, "divide by zero fault", "div.fault");
    }

    1000 / ( symVar - location);
#pragma clang optimize on
}
