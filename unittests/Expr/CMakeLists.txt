add_klee_unit_test(ExprTest
  ExprTest.cpp
  RangesTest.cpp)
target_link_libraries(ExprTest PRIVATE kleaverExpr)
