
#COMMON_FLAGS:=-external-calls=all -solver-backend=z3 -use-cex-cache=false -only-output-states-covering-new
COMMON_FLAGS:=-external-calls=all -solver-backend=z3 -only-output-states-covering-new
PREFIX:=
all: $(PREFIX)klee-out-int $(PREFIX)klee-out-normal
$(PREFIX)klee-out-int: srt.bc klee
	rm -rf $@
	./klee $(COMMON_FLAGS) -use-int-solver -debug-z3int-dump-queries=/tmp/z3int.smt2 -output-dir=$(PREFIX)$@ srt.bc

$(PREFIX)klee-out-normal: srt.bc
	rm -rf $@
	../../build7.0/bin/klee $(COMMON_FLAGS) -debug-z3-dump-queries=/tmp/z3.smt2 -output-dir=$(PREFIX)$@ srt.bc

$(PREFIX)klee-out-crosscheck: srt.bc klee
	rm -rf $@
	./klee $(COMMON_FLAGS) -use-int-solver -debug-crosscheck-z3int-solver=stp -debug-z3int-dump-queries=/tmp/z3int.smt2 -output-dir=$(PREFIX)$@ srt.bc
CU_FLAGS:=--disable-inlining --libc=uclibc --posix-runtime --switch-type=internal --use-batching-search --batch-instructions=10000 --simplify-sym-indices
$(PREFIX)klee-out-%: %.bc klee
	rm -rf $@
#	./klee $(COMMON_FLAGS) $(CU_FLAGS) -use-int-solver -output-dir=$@ $< -sym-args 0 1 10 --sym-args 0 2 2 --sym-files 1 8 --sym-stdin 8 --sym-stdout
	./klee $(COMMON_FLAGS) $(CU_FLAGS) -use-int-solver -debug-crosscheck-z3int-solver=stp -output-dir=$@ $< -sym-args 0 1 10 --sym-args 0 2 2 --sym-files 1 8 --sym-stdin 8 --sym-stdout
	#./klee $(COMMON_FLAGS) $(CU_FLAGS) -use-int-solver -debug-crosscheck-z3int-solver=stp -debug-assignment-validating-solver -output-dir=$@ $< -sym-args 0 1 10 --sym-args 0 2 2 --sym-files 1 8 --sym-stdin 8 --sym-stdout


srt.bc: ../../examples/sort/sort.c
	/data/dependencies/llvm-7.0/build/bin/clang -g -c -emit-llvm -I../../include -o $@ $<
