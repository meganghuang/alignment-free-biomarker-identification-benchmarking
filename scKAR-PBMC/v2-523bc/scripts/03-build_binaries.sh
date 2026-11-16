#!/bin/bash
# 03-build_binaries.sh — compile scKAR's C++ components.
# The Nov 12 run failed here: the build invoked a toolchain that could not
# resolve -lstdc++, so create_adj_list and ref_filter were never produced and
# every downstream stage silently emitted empty output. Build with g++ (which
# links libstdc++ itself) and fail loudly if a binary is missing.
set -uo pipefail
ROOT=/u/scratch/m/mhuang/claude/compare3/sckar_percell
SRC=$ROOT/repo/src
echo "=== build on $(hostname) $(date) ==="
g++ --version | head -1

build () {  # <out> <src...>
  local out="$1"; shift
  echo "--- $out ---"
  g++ -std=c++17 -O2 -o "$out" "$@" 2>&1 | tail -5
  if [ -x "$out" ]; then echo "OK $out"; else echo "FAILED $out"; return 1; fi
}

cd $SRC/kmer-filtering
build ./create_adj_list create_adj_list.cpp || exit 1
build ./ref_filter      ref_filter.cpp      || exit 1

cd $SRC/f-test
for f in adj_to_mat.cpp bh_correction.cpp f_oneway_test.cpp; do
  [ -f "$f" ] && build "./${f%.cpp}" "$f"
done
if [ -f f_test.cpp ]; then build ./f_test.out f_test.cpp; fi

cd $SRC/clustering
[ -f merge_adj_for_clustering.cpp ] && build ./merge_adj_for_clustering merge_adj_for_clustering.cpp

echo "=== built binaries ==="
find $SRC -maxdepth 2 -type f -perm -u+x ! -name "*.cpp" ! -name "*.sh" -printf "%p %s\n"
