// RUN: %llvmgcc %s -emit-llvm -g %O0opt -c -o %t.bc
// RUN: rm -rf %t.klee-out
// RUN: %klee --output-dir=%t.klee-out -exit-on-error-type=All %t.bc

#include <assert.h>
#include <klee/klee.h>

union U0 {
	signed f3 :18;
};

static union U0 u = { 0UL };

int main(int argc, char **argv) {
  u.f3 = 534;

  klee_assert(u.f3 == 534);

  return 0;
}

