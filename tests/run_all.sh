#!/usr/bin/env bash
# Unified test runner for VegasAfterglow.
#
#   tests/run_all.sh                 # everything: C++ + pytest + full validation + test-report.html
#   tests/run_all.sh --quick         # skip the validation suite (tiers 1+2 only, seconds)
#   tests/run_all.sh --build         # force-rebuild C++ tests and reinstall the module first
#   tests/run_all.sh --no-build      # skip the automatic staleness check
#   tests/run_all.sh --cpp           # C++ only
#   tests/run_all.sh --py            # pytest only
#   tests/run_all.sh --physics       # pytest -m "physics or golden" only
#   tests/run_all.sh --validation    # validation suite only (slow)
#   tests/run_all.sh --no-html       # skip the self-contained test-report.html
#   tests/run_all.sh --cov           # also write coverage summary + htmlcov/
#
# Anything the selected tiers need that is older than the sources (the C++ test
# binary, the installed Python module) is rebuilt automatically; testing a stale
# build silently is worse than the rebuild wait.
#
# Exit code is nonzero if any selected tier fails.
set -u

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

BOLD=$(tput bold 2>/dev/null || true)
RED=$(tput setaf 1 2>/dev/null || true)
GREEN=$(tput setaf 2 2>/dev/null || true)
CYAN=$(tput setaf 6 2>/dev/null || true)
RESET=$(tput sgr0 2>/dev/null || true)

RUN_CPP=1 RUN_PY=1 RUN_VALIDATION=1 DO_BUILD=0 NO_BUILD=0 DO_HTML=1 DO_COV=0 PY_MARK=""
for arg in "$@"; do
    case "$arg" in
        --quick) RUN_VALIDATION=0 ;;
        --cpp) RUN_PY=0 RUN_VALIDATION=0 ;;
        --py) RUN_CPP=0 RUN_VALIDATION=0 ;;
        --physics) RUN_CPP=0 RUN_VALIDATION=0 PY_MARK="physics or golden" ;;
        --validation) RUN_CPP=0 RUN_PY=0 RUN_VALIDATION=1 ;;
        --build) DO_BUILD=1 ;;
        --no-build) NO_BUILD=1 ;;
        --html) DO_HTML=1 ;;
        --no-html) DO_HTML=0 ;;
        --cov) DO_COV=1 ;;
        *) echo "unknown option: $arg"; exit 2 ;;
    esac
done

declare -a NAMES RESULTS TIMES
overall=0

run_tier() { # name, command...
    local name=$1; shift
    echo ""
    echo "${BOLD}${CYAN}=== $name ===${RESET}"
    local start=$SECONDS
    "$@"
    local rc=$?
    local dt=$((SECONDS - start))
    NAMES+=("$name"); TIMES+=("${dt}s")
    if [ $rc -eq 0 ]; then
        RESULTS+=("${GREEN}PASS${RESET}")
    else
        RESULTS+=("${RED}FAIL${RESET}"); overall=1
    fi
}

# A target is stale when any C++ source, build script, or packaging file is
# newer than it. `find -newer` prints the first offender (empty = fresh).
newer_than() { # ref-file; stdout: first source newer than it
    find src pybind tests/cpp CMakeLists.txt pyproject.toml \
        \( -name '*.cpp' -o -name '*.h' -o -name '*.hpp' -o -name '*.tpp' \
           -o -name 'CMakeLists.txt' -o -name 'pyproject.toml' \) \
        -newer "$1" -print 2>/dev/null | head -1
}

NEED_CPP_BUILD=$DO_BUILD
NEED_PY_BUILD=$DO_BUILD
if [ $NO_BUILD -eq 0 ] && [ $DO_BUILD -eq 0 ]; then
    if [ $RUN_CPP -eq 1 ]; then
        if [ ! -x build/afterglow_tests ]; then
            echo "build/afterglow_tests missing — will build"
            NEED_CPP_BUILD=1
        else
            OFFENDER=$(newer_than build/afterglow_tests)
            [ -n "$OFFENDER" ] && {
                echo "build/afterglow_tests is stale ($OFFENDER is newer) — rebuilding"
                NEED_CPP_BUILD=1
            }
        fi
    fi
    if [ $RUN_PY -eq 1 ] || [ $RUN_VALIDATION -eq 1 ]; then
        MODULE_SO=$(python -c "import VegasAfterglow.VegasAfterglowC as m; print(m.__file__)" 2>/dev/null)
        if [ -z "$MODULE_SO" ]; then
            echo "VegasAfterglow module not importable — will install"
            NEED_PY_BUILD=1
        else
            OFFENDER=$(newer_than "$MODULE_SO")
            [ -n "$OFFENDER" ] && {
                echo "installed module is stale ($OFFENDER is newer) — reinstalling"
                NEED_PY_BUILD=1
            }
        fi
    fi
fi

if [ $NEED_CPP_BUILD -eq 1 ] && [ $RUN_CPP -eq 1 ]; then
    run_tier "build: C++ tests" bash -c \
        "cmake -B build -DAFTERGLOW_TESTS=ON >/dev/null && cmake --build build -j8 | tail -1"
fi
if [ $NEED_PY_BUILD -eq 1 ] && { [ $RUN_PY -eq 1 ] || [ $RUN_VALIDATION -eq 1 ]; }; then
    run_tier "build: python module (pip install -e .)" pip install -e . -q
fi

if [ $RUN_CPP -eq 1 ]; then
    if [ ! -x build/afterglow_tests ]; then
        echo "${RED}build/afterglow_tests not found — run with --build first${RESET}"
        overall=1
    else
        CPP_ARGS=(--log_level=error)
        [ $DO_HTML -eq 1 ] && CPP_ARGS=("--logger=HRF,error:JUNIT,all,cpp-junit.xml")
        run_tier "C++ unit tests (Boost.Test)" ./build/afterglow_tests "${CPP_ARGS[@]}"
    fi
fi

if [ $RUN_PY -eq 1 ]; then
    echo "using python: $(command -v python)"
    python -c "import pytest, VegasAfterglow" 2>/dev/null || {
        echo "${RED}python env not ready:${RESET} needs pytest + VegasAfterglow installed"
        echo "  (activate the right environment, then: pip install -e '.[test]')"
        exit 2
    }
    PYTEST_ARGS=(tests/ -q)
    [ -n "$PY_MARK" ] && PYTEST_ARGS+=(-m "$PY_MARK")
    [ $DO_HTML -eq 1 ] && PYTEST_ARGS+=(--junitxml=py-junit.xml)
    if [ $DO_COV -eq 1 ]; then
        python -c "import pytest_cov" 2>/dev/null \
            || { echo "${RED}--cov needs pytest-cov:${RESET} pip install pytest-cov"; exit 2; }
        PYTEST_ARGS+=(--cov=VegasAfterglow --cov-report=term --cov-report=html)
    fi
    run_tier "Python tests (pytest)" python -m pytest "${PYTEST_ARGS[@]}"
fi

if [ $RUN_VALIDATION -eq 1 ]; then
    run_tier "Full validation suite (slow)" python tests/validation/run_validation.py --all
fi

if [ $DO_HTML -eq 1 ]; then
    REPORT_XML=()
    REPORT_LABELS=()
    [ -f cpp-junit.xml ] && REPORT_XML+=(cpp-junit.xml) \
        && REPORT_LABELS+=(--label "cpp-junit.xml=C++ (Boost.Test)")
    [ -f py-junit.xml ] && REPORT_XML+=(py-junit.xml) \
        && REPORT_LABELS+=(--label "py-junit.xml=Python (pytest)")
    VAL_JSON=tests/validation/regression/results/regression_results.json
    BENCH_JSON=tests/validation/benchmark/results/benchmark_history.json
    VAL_ARG=()
    [ -f "$VAL_JSON" ] && VAL_ARG+=(--validation "$VAL_JSON")
    [ -f "$BENCH_JSON" ] && VAL_ARG+=(--benchmark "$BENCH_JSON")
    if [ ${#REPORT_XML[@]} -gt 0 ]; then
        python tests/report.py "${REPORT_XML[@]}" "${REPORT_LABELS[@]}" \
            ${VAL_ARG[@]+"${VAL_ARG[@]}"} -o test-report.html || true
    fi
fi

echo ""
echo "${BOLD}=== Summary ===${RESET}"
[ ${#NAMES[@]} -eq 0 ] && { echo "  (no tiers ran)"; exit $overall; }
for i in "${!NAMES[@]}"; do
    printf "  %-42s %b  (%s)\n" "${NAMES[$i]}" "${RESULTS[$i]}" "${TIMES[$i]}"
done
[ $DO_HTML -eq 1 ] && echo "  report: test-report.html (unified: tests + validation + performance)"
[ $DO_COV -eq 1 ] && echo "  coverage: htmlcov/index.html"
exit $overall
