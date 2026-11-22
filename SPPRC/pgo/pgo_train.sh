#!/bin/bash
# PGO training for libSPPRC (see SPPRC_PGO_GEN / SPPRC_PGO_USE in CMakeLists.txt).
#
#  1. builds the instrumented library in a temporary folder (cmake options copied from
#     build/CMakeCache.txt: Eigen, Boost, OpenMP);
#  2. runs the root node of the training instances with Coluna (pgo/train.jl);
#  3. merges the profiles into pgo/default.profdata.
#
# Then build the optimized library:
#   cmake -S . -B build -DSPPRC_PGO_USE=pgo/default.profdata && ninja -C build
#
# Training: roots of E-n51-k5, E-n76-k7, E-n101-k14 (E-n101-k8 is left out to test the
# profile on an instance it did not see). Takes about 2 minutes.
# A short tree can be added to TRAIN (third field = seconds, e.g. "E-n76-k8 736 60"),
# but it did not help (2026-11-22, 4-min trees E-n76-k7 / E-n101-k8: same labeling
# times as with the roots only).
#
# The instrumented build is made in $TMPDIR, not inside the project: inside Claude
# Code's sandbox, CMake cannot identify the compiler in a new build folder in this
# project (the compiler is then "unknown" and the -std=gnu++23 flag is lost). In a
# normal terminal both places work.
#
# Environment (optional): JULIA, SYSIMAGE (default ~/coluna_sys.dylib, "" = none),
# CMAKE, NINJA, PROFDATA (llvm-profdata of the same clang as the build).

set -euo pipefail

ROOT=$(cd "$(dirname "$0")/.." && pwd)
OUT="$ROOT/pgo/default.profdata"
CACHE="$ROOT/build/CMakeCache.txt"

JULIA=${JULIA:-julia}
SYSIMAGE=${SYSIMAGE-$HOME/coluna_sys.dylib}
CMAKE=${CMAKE:-$(grep -m1 '^CMAKE_COMMAND:' "$CACHE" 2>/dev/null | cut -d= -f2- || true)}
CMAKE=${CMAKE:-cmake}
NINJA=${NINJA:-$(grep -m1 '^CMAKE_MAKE_PROGRAM:' "$CACHE" 2>/dev/null | cut -d= -f2- || true)}
NINJA=${NINJA:-ninja}
PROFDATA=${PROFDATA:-$(xcrun -f llvm-profdata 2>/dev/null || command -v llvm-profdata)}

# instance, primal bound (opt + 1), [tree seconds] (none = root only)
TRAIN=("E-n51-k5 522" "E-n76-k7 683" "E-n101-k14 1068")

# Work folder: build/ (instrumented library, loaded as ../build/libSPPRC.dylib from
# run/, see Coluna/instanceCVRP.jl), run/, prof/ (.profraw files), logs.
WORK=$(mktemp -d "${TMPDIR:-/tmp}/spprc_pgo.XXXXXX")
GEN="$WORK/build"
mkdir -p "$WORK/run" "$WORK/prof"

# 1. instrumented build, with the same Eigen/Boost/OpenMP options as build/
CMAKE_ARGS=(-G Ninja "-DCMAKE_MAKE_PROGRAM=$NINJA" -DCMAKE_BUILD_TYPE=Release -DSPPRC_PGO_GEN=ON)
if [ -f "$CACHE" ]; then
    while IFS= read -r line; do
        CMAKE_ARGS+=("-D${line%%:*}=${line#*=}")
    done < <(grep -E '^(Eigen3_DIR|Boost_DIR|OpenMP_(CXX|C)_(FLAGS|LIB_NAMES|INCLUDE_DIR)|OpenMP_libomp_LIBRARY):' "$CACHE")
fi
echo "== instrumented build in $GEN"
"$CMAKE" -S "$ROOT" -B "$GEN" "${CMAKE_ARGS[@]}" > "$WORK/cmake.log" 2>&1 \
    || { cat "$WORK/cmake.log"; exit 1; }
"$NINJA" -C "$GEN" > "$WORK/ninja.log" 2>&1 || { tail -30 "$WORK/ninja.log"; exit 1; }

# 2. training runs
JFLAGS=()
[ -n "$SYSIMAGE" ] && JFLAGS=(-J "$SYSIMAGE")
for t in "${TRAIN[@]}"; do
    set -- $t
    if [ -n "${3:-}" ]; then echo "== training: $1 tree, $3 s"; else echo "== training: $1 root"; fi
    ( cd "$WORK/run" && LLVM_PROFILE_FILE="$WORK/prof/$1-%p.profraw" \
        "$JULIA" ${JFLAGS[@]+"${JFLAGS[@]}"} "$ROOT/pgo/train.jl" "$1" "$2" ${3:+"$3"} > "$WORK/$1.log" 2>&1 ) \
        || { tail -30 "$WORK/$1.log"; exit 1; }
    grep -E "Dual bound" "$WORK/$1.log" || true
done

# 3. merge
"$PROFDATA" merge -o "$OUT" "$WORK"/prof/*.profraw
echo "== profile written: $OUT"
echo "   build:  cmake -S \"$ROOT\" -B \"$ROOT/build\" -DSPPRC_PGO_USE=pgo/default.profdata && ninja -C \"$ROOT/build\""
rm -rf "$WORK"
