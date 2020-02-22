#!/usr/bin/env python3
import sys
import os
import subprocess
import tempfile

GEN_BOUT="/data/klee/build7.0/bin/gen-bout"
KLEE="/data/klee/build7.0/bin/klee"
#KLEE="klee"
#GEN_BOUT="gen-bout"

def split_args():
  prog = None
  prog_args = []
  klee_args = []
  is_progargs = False
  for a in sys.argv[1:]:
      if is_progargs:
          prog_args += [a]
      elif a.startswith("-"):
          klee_args += [a]
      else:
          prog = a
          is_progargs = True
  return klee_args, prog, prog_args

def maybe_file_size(name):
  try:
    with open(name, 'r') as f:
      f.seek(0, os.SEEK_END)
      return f.tell()
  except:
    return None

def prog_args_to_posix(prog_args):
  posix_args = []
  sym_file = 'A'
  sym_file_sizes = [] 
  gen_out_args = []
  for parg in prog_args:
      file_size = maybe_file_size(parg)
      if file_size is None:
          posix_args += ['--sym-arg', str(len(parg))]
          gen_out_args += [parg]
      else:
          sym_file_sizes += [file_size]
          posix_args += [sym_file]
          sym_file = chr(ord(sym_file) + 1)
          gen_out_args += ['--sym-file', parg]

  if ord(sym_file) - ord('A') > 0:
      posix_args += ['--sym-files', str(ord(sym_file) - ord('A')), str(max(sym_file_sizes))]
  return posix_args, gen_out_args

def create_ktest_file(gen_out_args, tmpdir):
  out_file=tmpdir + "/test.ktest"
  subprocess.run([GEN_BOUT, "--bout-file", out_file] + gen_out_args, check=True)
  return out_file

def main():
  tmpdir = tempfile.TemporaryDirectory()
  klee_args, prog, prog_args = split_args()
  posix_args, gen_out_args = prog_args_to_posix(prog_args)
  ktest_file = create_ktest_file(gen_out_args,tmpdir.name)
  klee_args += ["-seed-file=" + ktest_file]
#  print(subprocess.run(["/data/klee/build7.0/bin/ktest-tool", ktest_file], stdout=sys.stdout))
  subprocess.run([KLEE] + klee_args + [prog] + posix_args, stdout=sys.stdout)
#  print(subprocess.run(["/data/klee/build7.0/bin/ktest-tool", ktest_file], stdout=sys.stdout))
  print(posix_args)
  print(ktest_file)
  print(prog_args)




main()


