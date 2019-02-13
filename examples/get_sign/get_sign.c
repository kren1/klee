/*
 * First KLEE tutorial: testing a small function
 */

#include <klee/klee.h>

int get_sign(int x) {
  if (x == -1)
     return 0;
  
  if (x < -2)
     return -1;
  else 
     return 1;
} 

int main() {
  int a;
  klee_make_symbolic(&a, sizeof(a), "a");
  return get_sign(a);
} 
