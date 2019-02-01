#===------------------------------------------------------------------------===#
#
#                     The KLEE Symbolic Virtual Machine
#
# This file is distributed under the University of Illinois Open Source
# License. See LICENSE.TXT for details.
#
#===------------------------------------------------------------------------===#

## Build 64bit library
add_library(kleeRuntest64 SHARED
  intrinsics.c
  KTest.cpp
)
# Increment version appropriately if ABI/API changes, more details:
# http://tldp.org/HOWTO/Program-Library-HOWTO/shared-libraries.html#AEN135
set(KLEE_RUNTEST_VERSION 1.0)
set_target_properties(kleeRuntest64
  PROPERTIES
    VERSION ${KLEE_RUNTEST_VERSION}
    SOVERSION ${KLEE_RUNTEST_VERSION}
)

install(TARGETS kleeRuntest64
  DESTINATION "${CMAKE_INSTALL_FULL_LIBDIR}")

add_library(kleeRuntest32 SHARED
  intrinsics.c
  KTest.cpp
)

## Build 64bit library
target_compile_options(kleeRuntest32
	PRIVATE -m32
	)

target_link_libraries(kleeRuntest32
	-m32
	)

# Increment version appropriately if ABI/API changes, more details:
# http://tldp.org/HOWTO/Program-Library-HOWTO/shared-libraries.html#AEN135
set(KLEE_RUNTEST_VERSION 1.0)
set_target_properties(kleeRuntest32
  PROPERTIES
    VERSION ${KLEE_RUNTEST_VERSION}
    SOVERSION ${KLEE_RUNTEST_VERSION}
)

install(TARGETS kleeRuntest32
  DESTINATION "${CMAKE_INSTALL_FULL_LIBDIR}")
