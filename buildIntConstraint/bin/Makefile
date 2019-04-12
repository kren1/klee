
#COMMON_FLAGS:=-external-calls=all -solver-backend=z3 -use-cex-cache=false -only-output-states-covering-new
COMMON_FLAGS:=-external-calls=all -solver-backend=z3 -only-output-states-covering-new
PREFIX:=

inc: $(foreach var, $(shell seq 1 8), srt-$(var))

srt-%: $(PREFIX)klee-out-int-% $(PREFIX)klee-out-noint-%
	true

$(PREFIX)klee-out-int-%: srt.%.bc klee
	rm -rf $@
	./klee $(COMMON_FLAGS) -use-int-solver -output-dir=$@ $<

$(PREFIX)klee-out-noint-%: srt.%.bc klee
	rm -rf $@
	./klee $(COMMON_FLAGS) -output-dir=$@ $<


$(PREFIX)klee-out-count: srt.bc klee
	rm -rf $@
	./klee $(COMMON_FLAGS) -count-constructs -use-query-log=solver:kquery -use-int-solver -output-dir=$@ srt.bc



$(PREFIX)klee-out-normal: srt.bc
	rm -rf $@
	../../build7.0/bin/klee $(COMMON_FLAGS) -output-dir=$@ srt.bc

$(PREFIX)klee-out-crosscheck: srt.bc klee
	rm -rf $@
	./klee $(COMMON_FLAGS) -use-int-solver -debug-crosscheck-z3int-solver=stp -debug-z3int-dump-queries=/tmp/z3int.smt2 -output-dir=$(PREFIX)$@ srt.bc
CU_FLAGS:=--disable-inlining --libc=uclibc --posix-runtime --switch-type=internal --use-batching-search --batch-instructions=10000 --simplify-sym-indices
$(PREFIX)klee-out-%: %.bc klee
	rm -rf $@
#	./klee $(COMMON_FLAGS) $(CU_FLAGS) -use-int-solver -output-dir=$@ $< -sym-args 0 1 10 --sym-args 0 2 2 --sym-files 1 8 --sym-stdin 8 --sym-stdout
	./klee $(COMMON_FLAGS) $(CU_FLAGS) -use-query-log=all:kquery -use-int-solver -debug-crosscheck-z3int-solver=stp -output-dir=$@ $< -sym-args 0 1 10 --sym-args 0 2 2 --sym-files 1 8 --sym-stdin 8 --sym-stdout
	#./klee $(COMMON_FLAGS) $(CU_FLAGS) -use-int-solver -debug-crosscheck-z3int-solver=stp -debug-assignment-validating-solver -output-dir=$@ $< -sym-args 0 1 10 --sym-args 0 2 2 --sym-files 1 8 --sym-stdin 8 --sym-stdout


srt.%.bc: ../../examples/sort/sort.c
	/data/dependencies/llvm-7.0/build/bin/clang -DN=$* -g -c -emit-llvm -I../../include -o $@ $<
