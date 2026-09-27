#!/usr/bin/env bash
# fTetWild CMake configuration script.
#
# Usage: bash configure.sh [options] [-- extra cmake args]
#
# All options have sensible defaults; run with --help for the full list.

set -euo pipefail

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
RED='\033[0;31m'; YELLOW='\033[1;33m'; GREEN='\033[0;32m'; CYAN='\033[0;36m'
BOLD='\033[1m'; RESET='\033[0m'

info()  { echo -e "${GREEN}[configure]${RESET} $*"; }
warn()  { echo -e "${YELLOW}[configure] WARNING:${RESET} $*" >&2; }
error() { echo -e "${RED}[configure] ERROR:${RESET} $*" >&2; exit 1; }

usage() {
cat <<EOF
${BOLD}Usage:${RESET} bash configure.sh [options] [-- cmake-extra-args...]

${BOLD}Build options:${RESET}
  --build-type TYPE         CMake build type: Release (default), Debug,
                            RelWithDebInfo, MinSizeRel
  --generator GEN           CMake generator: Ninja (default if available),
                            "Unix Makefiles"
  --install-prefix PATH     CMAKE_INSTALL_PREFIX  [default: /usr/local]
  --jobs N                  Hint printed for the subsequent build step

${BOLD}Compiler options:${RESET}
  --cxx COMPILER            C++ compiler (e.g. g++, clang++)
  --cc  COMPILER            C compiler   (e.g. gcc,  clang)
  --cxx-flags FLAGS         Extra CMAKE_CXX_FLAGS (quoted)

${BOLD}fTetWild feature flags:${RESET}
  --no-tbb                  Disable TBB parallelism      [default: ON]
  --use-float               Use float instead of double  [default: OFF]
  --exact-envelope          Enable FastEnvelope library  [default: OFF]
  --with-tetgen             Enable TetGen via libigl     [default: OFF]
  --with-gmp                Use GMP for exact rational AMIPS energy
                            stabilisation. Requires libgmp-dev.
                            Default: OFF (long double fallback is used)

${BOLD}Sanitizers (implies --build-type=RelWithDebInfo unless set):${RESET}
  --sanitize-address        Enable AddressSanitizer   (ASan)
  --sanitize-memory         Enable MemorySanitizer    (MSan) — clang only
  --sanitize-thread         Enable ThreadSanitizer    (TSan)
  --sanitize-undefined      Enable UndefinedBehaviorSanitizer (UBSan)
  Note: ASan and TSan/MSan are mutually exclusive.

${BOLD}Test / output options:${RESET}
  --no-tests                Disable unit test target     [default: ON]
  --verbose                 CMAKE_VERBOSE_MAKEFILE=ON
  --no-compile-commands     Disable compile_commands.json generation

${BOLD}Dependency / path options:${RESET}
  --gmp-inc PATH            GMP include directory (if not auto-detected)
  --gmp-lib PATH            GMP library directory (if not auto-detected)
  --third-party-dir DIR     Alternate directory for 3rd-party sources
  --offline                 FETCHCONTENT_UPDATES_DISCONNECTED=ON
                            (skip network checks for already-fetched deps)

${BOLD}Misc:${RESET}
  --dry-run                 Print the cmake command without running it
  -h, --help                Show this help message

${BOLD}Examples:${RESET}
  bash configure.sh
  bash configure.sh --build-type Debug --sanitize-address
  bash configure.sh --build-type RelWithDebInfo --sanitize-undefined --cxx clang++
  bash configure.sh --exact-envelope --no-tbb
  bash configure.sh --with-gmp                   # enable exact GMP rational arithmetic
  bash configure.sh --offline -- -DCMAKE_VERBOSE_MAKEFILE=ON
EOF
}

# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------
BUILD_TYPE="Release"
GENERATOR=""                   # auto-detected below
INSTALL_PREFIX=""
JOBS=""

CXX_COMPILER=""
C_COMPILER=""
CXX_FLAGS=""

ENABLE_TBB=ON
USE_FLOAT=OFF
EXACT_ENVELOPE=OFF
WITH_TETGEN=OFF
USE_GMP=OFF

WITH_SANITIZERS=OFF
SANITIZE_ADDRESS=OFF
SANITIZE_MEMORY=OFF
SANITIZE_THREAD=OFF
SANITIZE_UNDEFINED=OFF

BUILD_TESTING=ON
VERBOSE=OFF
EXPORT_COMPILE_COMMANDS=ON

GMP_INC=""
GMP_LIB=""
THIRD_PARTY_DIR=""
OFFLINE=OFF

DRY_RUN=OFF

EXTRA_ARGS=()

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------
while [[ $# -gt 0 ]]; do
    case "$1" in
        --build-type)        BUILD_TYPE="$2";        shift 2 ;;
        --build-type=*)      BUILD_TYPE="${1#*=}";   shift   ;;
        --generator)         GENERATOR="$2";         shift 2 ;;
        --generator=*)       GENERATOR="${1#*=}";    shift   ;;
        --install-prefix)    INSTALL_PREFIX="$2";    shift 2 ;;
        --install-prefix=*)  INSTALL_PREFIX="${1#*=}"; shift ;;
        --jobs)              JOBS="$2";              shift 2 ;;
        --jobs=*)            JOBS="${1#*=}";         shift   ;;

        --cxx)               CXX_COMPILER="$2";      shift 2 ;;
        --cxx=*)             CXX_COMPILER="${1#*=}"; shift   ;;
        --cc)                C_COMPILER="$2";        shift 2 ;;
        --cc=*)              C_COMPILER="${1#*=}";   shift   ;;
        --cxx-flags)         CXX_FLAGS="$2";         shift 2 ;;
        --cxx-flags=*)       CXX_FLAGS="${1#*=}";    shift   ;;

        --no-tbb)            ENABLE_TBB=OFF;         shift   ;;
        --use-float)         USE_FLOAT=ON;           shift   ;;
        --exact-envelope)    EXACT_ENVELOPE=ON;      shift   ;;
        --with-tetgen)       WITH_TETGEN=ON;         shift   ;;
        --with-gmp)          USE_GMP=ON;             shift   ;;

        --sanitize-address)  SANITIZE_ADDRESS=ON;   WITH_SANITIZERS=ON; shift ;;
        --sanitize-memory)   SANITIZE_MEMORY=ON;    WITH_SANITIZERS=ON; shift ;;
        --sanitize-thread)   SANITIZE_THREAD=ON;    WITH_SANITIZERS=ON; shift ;;
        --sanitize-undefined) SANITIZE_UNDEFINED=ON; WITH_SANITIZERS=ON; shift ;;

        --no-tests)          BUILD_TESTING=OFF;      shift   ;;
        --verbose)           VERBOSE=ON;             shift   ;;
        --no-compile-commands) EXPORT_COMPILE_COMMANDS=OFF; shift ;;

        --gmp-inc)           GMP_INC="$2";           shift 2 ;;
        --gmp-inc=*)         GMP_INC="${1#*=}";      shift   ;;
        --gmp-lib)           GMP_LIB="$2";           shift 2 ;;
        --gmp-lib=*)         GMP_LIB="${1#*=}";      shift   ;;
        --third-party-dir)   THIRD_PARTY_DIR="$2";   shift 2 ;;
        --third-party-dir=*) THIRD_PARTY_DIR="${1#*=}"; shift ;;
        --offline)           OFFLINE=ON;             shift   ;;

        --dry-run)           DRY_RUN=ON;             shift   ;;
        -h|--help)           usage; exit 0 ;;

        --)  shift; EXTRA_ARGS+=("$@"); break ;;
        -D*) EXTRA_ARGS+=("$1"); shift ;;           # raw -D pass-through
        *)   error "Unknown option: $1\nRun with --help for usage." ;;
    esac
done

# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------
valid_build_types="Release Debug RelWithDebInfo MinSizeRel"
if [[ ! " $valid_build_types " =~ " $BUILD_TYPE " ]]; then
    error "Invalid --build-type '$BUILD_TYPE'. Valid: $valid_build_types"
fi

if [[ "$SANITIZE_ADDRESS" == "ON" && ("$SANITIZE_THREAD" == "ON" || "$SANITIZE_MEMORY" == "ON") ]]; then
    error "ASan is incompatible with TSan/MSan. Choose one set."
fi

if [[ "$WITH_SANITIZERS" == "ON" && "$BUILD_TYPE" == "Release" ]]; then
    warn "Sanitizers are enabled but build type is Release."
    warn "Switching to RelWithDebInfo for useful stack traces."
    BUILD_TYPE="RelWithDebInfo"
fi

if [[ "$SANITIZE_MEMORY" == "ON" ]]; then
    if ! "${CXX_COMPILER:-c++}" --version 2>&1 | grep -qi clang; then
        warn "MSan is only reliable with Clang. Consider --cxx clang++."
    fi
fi

if [[ "$USE_FLOAT" == "ON" && "$EXACT_ENVELOPE" == "ON" ]]; then
    warn "--use-float with --exact-envelope has not been tested; proceed with caution."
fi

# ---------------------------------------------------------------------------
# Auto-detect generator
# ---------------------------------------------------------------------------
if [[ -z "$GENERATOR" ]]; then
    if command -v ninja &>/dev/null; then
        GENERATOR="Ninja"
    else
        GENERATOR="Unix Makefiles"
    fi
fi

# ---------------------------------------------------------------------------
# Locate source directory (one level above this script)
# ---------------------------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SOURCE_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

# ---------------------------------------------------------------------------
# Assemble cmake arguments
# ---------------------------------------------------------------------------
CMAKE_ARGS=(
    -G "${GENERATOR}"
    -DCMAKE_BUILD_TYPE="${BUILD_TYPE}"
    -DCMAKE_EXPORT_COMPILE_COMMANDS="${EXPORT_COMPILE_COMMANDS}"
    -DBUILD_TESTING="${BUILD_TESTING}"
    -DCMAKE_VERBOSE_MAKEFILE="${VERBOSE}"

    # fTetWild feature flags
    -DFLOAT_TETWILD_ENABLE_TBB="${ENABLE_TBB}"
    -DFLOAT_TETWILD_USE_FLOAT="${USE_FLOAT}"
    -DFLOAT_TETWILD_WITH_EXACT_ENVELOPE="${EXACT_ENVELOPE}"
    -DFLOAT_TETWILD_WITH_SANITIZERS="${WITH_SANITIZERS}"
    -DFLOAT_TETWILD_USE_GMP="${USE_GMP}"

    # Sanitizers
    -DSANITIZE_ADDRESS="${SANITIZE_ADDRESS}"
    -DSANITIZE_MEMORY="${SANITIZE_MEMORY}"
    -DSANITIZE_THREAD="${SANITIZE_THREAD}"
    -DSANITIZE_UNDEFINED="${SANITIZE_UNDEFINED}"

    # libigl extras
    -DLIBIGL_WITH_TETGEN="${WITH_TETGEN}"
    -DLIBIGL_WITH_PREDICATES=ON

    # FetchContent
    -DFETCHCONTENT_UPDATES_DISCONNECTED="${OFFLINE}"
)

[[ -n "$INSTALL_PREFIX"  ]] && CMAKE_ARGS+=(-DCMAKE_INSTALL_PREFIX="${INSTALL_PREFIX}")
[[ -n "$CXX_COMPILER"   ]] && CMAKE_ARGS+=(-DCMAKE_CXX_COMPILER="${CXX_COMPILER}")
[[ -n "$C_COMPILER"     ]] && CMAKE_ARGS+=(-DCMAKE_C_COMPILER="${C_COMPILER}")
[[ -n "$CXX_FLAGS"      ]] && CMAKE_ARGS+=(-DCMAKE_CXX_FLAGS="${CXX_FLAGS}")
[[ -n "$GMP_INC"        ]] && CMAKE_ARGS+=(-DGMP_INC="${GMP_INC}")
[[ -n "$GMP_LIB"        ]] && CMAKE_ARGS+=(-DGMP_LIB="${GMP_LIB}")
[[ -n "$THIRD_PARTY_DIR" ]] && CMAKE_ARGS+=(-DINPUT_THIRD_PARTY_DIR="${THIRD_PARTY_DIR}")

# User pass-through
CMAKE_ARGS+=("${EXTRA_ARGS[@]}")

# ---------------------------------------------------------------------------
# Print configuration summary
# ---------------------------------------------------------------------------
echo ""
echo -e "${BOLD}${CYAN}━━━  fTetWild configuration  ━━━${RESET}"
printf "  %-28s %s\n" "Source dir:"        "${SOURCE_DIR}"
printf "  %-28s %s\n" "Build dir:"         "$(pwd)"
printf "  %-28s %s\n" "Generator:"         "${GENERATOR}"
printf "  %-28s %s\n" "Build type:"        "${BUILD_TYPE}"
[[ -n "$INSTALL_PREFIX" ]] && \
printf "  %-28s %s\n" "Install prefix:"    "${INSTALL_PREFIX}"
echo ""
printf "  %-28s %s\n" "TBB parallelism:"   "${ENABLE_TBB}"
printf "  %-28s %s\n" "Scalar type:"       "$( [[ $USE_FLOAT == ON ]] && echo float || echo double )"
printf "  %-28s %s\n" "Exact envelope:"    "${EXACT_ENVELOPE}"
printf "  %-28s %s\n" "TetGen:"            "${WITH_TETGEN}"
printf "  %-28s %s\n" "GMP rational:"      "$( [[ $USE_GMP == ON ]] && echo "ON (exact)" || echo "OFF (long double fallback)" )"
printf "  %-28s %s\n" "Build tests:"       "${BUILD_TESTING}"
printf "  %-28s %s\n" "compile_commands:"  "${EXPORT_COMPILE_COMMANDS}"
printf "  %-28s %s\n" "Verbose build:"     "${VERBOSE}"
printf "  %-28s %s\n" "Offline deps:"      "${OFFLINE}"
echo ""
if [[ "$WITH_SANITIZERS" == "ON" ]]; then
    SANS=()
    [[ $SANITIZE_ADDRESS   == ON ]] && SANS+=(ASan)
    [[ $SANITIZE_MEMORY    == ON ]] && SANS+=(MSan)
    [[ $SANITIZE_THREAD    == ON ]] && SANS+=(TSan)
    [[ $SANITIZE_UNDEFINED == ON ]] && SANS+=(UBSan)
    printf "  %-28s %s\n" "Sanitizers:" "${SANS[*]}"
    echo ""
fi
[[ -n "$CXX_COMPILER" ]] && printf "  %-28s %s\n" "C++ compiler:" "${CXX_COMPILER}"
[[ -n "$C_COMPILER"   ]] && printf "  %-28s %s\n" "C compiler:"   "${C_COMPILER}"
[[ -n "$CXX_FLAGS"    ]] && printf "  %-28s %s\n" "Extra CXX flags:" "${CXX_FLAGS}"
[[ -n "$GMP_INC"      ]] && printf "  %-28s %s\n" "GMP include:" "${GMP_INC}"
[[ -n "$GMP_LIB"      ]] && printf "  %-28s %s\n" "GMP lib dir:" "${GMP_LIB}"
[[ ${#EXTRA_ARGS[@]}  -gt 0 ]] && \
printf "  %-28s %s\n" "Extra cmake args:" "${EXTRA_ARGS[*]}"
echo ""

# ---------------------------------------------------------------------------
# Run (or print) cmake
# ---------------------------------------------------------------------------
CMD=(cmake "${SOURCE_DIR}" "${CMAKE_ARGS[@]}")

if [[ "$DRY_RUN" == "ON" ]]; then
    echo -e "${YELLOW}[dry-run]${RESET} Would execute:"
    echo "  ${CMD[*]}"
    echo ""
    exit 0
fi

info "Running cmake..."
echo ""
"${CMD[@]}"

echo ""
info "Configuration complete."
if [[ -n "$JOBS" ]]; then
    info "Build with:  cmake --build . --parallel ${JOBS}"
else
    info "Build with:  cmake --build . --parallel \$(nproc)"
fi
[[ "$BUILD_TESTING" == "ON" ]] && info "Test  with:  ctest --output-on-failure"
echo ""
