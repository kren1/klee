
#COMMON_FLAGS:=-external-calls=all -solver-backend=z3 -use-cex-cache=false -only-output-states-covering-new
COMMON_FLAGS:=-external-calls=all -solver-backend=z3 -only-output-states-covering-new
PREFIX:=
all: $(PREFIX)klee-out-int $(PREFIX)klee-out-normal
$(PREFIX)klee-out-int: srt.bc
	rm -rf $@
	./klee $(COMMON_FLAGS) -debug-z3int-dump-queries=/tmp/z3int.smt2 -output-dir=$(PREFIX)$@ srt.bc

$(PREFIX)klee-out-normal: srt.bc
	rm -rf $@
	../../build7.0/bin/klee $(COMMON_FLAGS) -debug-z3-dump-queries=/tmp/z3.smt2 -output-dir=$(PREFIX)$@ srt.bc

srt.bc: ../../examples/sort/sort.c
	/data/dependencies/llvm-7.0/build/bin/clang -g -c -emit-llvm -I../../include -o $@ $<
