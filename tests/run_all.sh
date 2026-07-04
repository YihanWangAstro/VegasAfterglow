#!/usr/bin/env bash
# Unified test runner for VegasAfterglow.
#
#   tests/run_all.sh                 # everything: C++ + pytest + full validation + test-report.html
#   tests/run_all.sh --quick         # skip the validation suite (tiers 1+2 only, seconds)
#   tests/run_all.sh --build         # force-reinstall the python module first
#   tests/run_all.sh --no-build      # skip the automatic build/freshness step
#   tests/run_all.sh --cpp           # C++ only
#   tests/run_all.sh --py            # pytest only
#   tests/run_all.sh --physics       # pytest -m "physics or golden" only
#   tests/run_all.sh --validation    # validation suite only (slow)
#   tests/run_all.sh --no-html       # skip the self-contained test-report.html
#   tests/run_all.sh --cov           # also write coverage summary + htmlcov/
#
# The build is kept fresh automatically — testing a stale build silently is
# worse than the rebuild wait. The C++ tier defers to cmake, whose dependency
# tracking makes a fresh tree a near-instant no-op; the installed Python module
# has no incremental mechanism, so it is reinstalled when it is missing, older
# than the sources, or built without profiling. The module is always built with
# AFTERGLOW_PROFILE=ON (like the deploy workflow) so the report's per-stage
# performance breakdown works.
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

if [ $RUN_CPP -eq 1 ] && [ $NO_BUILD -eq 0 ]; then
    # cmake's own dependency tracking is the freshness oracle: a fresh tree is a
    # near-instant no-op, and the heuristic alternatives (mtime scans) can only
    # drift from the real dependency graph.
    run_tier "build: C++ tests" bash -c \
        "cmake -B build -DAFTERGLOW_TESTS=ON >/dev/null && cmake --build build -j8 --target afterglow_tests | tail -1"
fi

if { [ $RUN_PY -eq 1 ] || [ $RUN_VALIDATION -eq 1 ]; } && [ $NO_BUILD -eq 0 ]; then
    NEED_PY_BUILD=$DO_BUILD
    if [ $NEED_PY_BUILD -eq 0 ]; then
        # One probe answers both freshness questions: where the installed module
        # lives, and whether it was built with profiling (probe the symbol the
        # report actually consumes).
        PROBE=$(python -c "import VegasAfterglow as va
print(va.VegasAfterglowC.__file__)
print(int(hasattr(va.Model, 'profile_data')))" 2>/dev/null)
        MODULE_SO=$(printf '%s\n' "$PROBE" | sed -n 1p)
        HAS_PROFILING=$(printf '%s\n' "$PROBE" | sed -n 2p)
        OFFENDER=""
        [ -n "$MODULE_SO" ] && OFFENDER=$(find src pybind external CMakeLists.txt pyproject.toml \
            \( -name '*.cpp' -o -name '*.h' -o -name '*.hpp' -o -name '*.tpp' \
               -o -name 'CMakeLists.txt' -o -name 'pyproject.toml' \) \
            -newer "$MODULE_SO" -print -quit)
        if [ -z "$MODULE_SO" ]; then
            echo "VegasAfterglow module not importable — will install"
            NEED_PY_BUILD=1
        elif [ -n "$OFFENDER" ]; then
            echo "installed module is stale ($OFFENDER is newer) — reinstalling"
            NEED_PY_BUILD=1
        elif [ "$HAS_PROFILING" != 1 ]; then
            echo "installed module lacks profiling (the report's per-stage timings need it) — reinstalling"
            NEED_PY_BUILD=1
        fi
    fi
    if [ $NEED_PY_BUILD -eq 1 ]; then
        run_tier "build: python module (pip install -e ., profiling on)" \
            pip install -e . -q --config-settings=cmake.define.AFTERGLOW_PROFILE=ON
    fi
fi

if [ $RUN_CPP -eq 1 ]; then
    if [ ! -x build/afterglow_tests ]; then
        echo "${RED}build/afterglow_tests not found — the build failed or --no-build was given${RESET}"
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
