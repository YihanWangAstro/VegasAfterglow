# === Platform detection ===
ifeq ($(OS),Windows_NT)
    PLATFORM := Windows
    EXE_EXT := .exe
else
    UNAME_S := $(shell uname -s)
    ifeq ($(UNAME_S),Darwin)
        PLATFORM := macOS
    else
        PLATFORM := Linux
    endif
    EXE_EXT :=
endif

# === Paths & sources ===
SRC_DIR := src
TEST_DIR := tests

# === Compiler & flags ===
CXX ?= g++
CXXFLAGS := -std=c++20 -O3 -march=native -flto -Iinclude -Isrc -Iexternal -w -g -DNDEBUG -DXTENSOR_DISABLE_EXCEPTIONS -DXTENSOR_USE_XSIMD
LDFLAGS := -lz

# Tests: each tests/foo/bar.cpp → tests/foo/bar(.exe)
TEST_SRCS := $(shell find $(TEST_DIR) -type f -name '*.cpp')
TEST_EXES := $(patsubst %.cpp,%$(EXE_EXT),$(TEST_SRCS))

# === Default target ===
.PHONY: all tests clean test

# Default to building tests
all: tests

# Individual targets
tests: $(TEST_EXES)

# Test executables - compile directly with source files
%$(EXE_EXT): %.cpp
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $< $(shell find $(SRC_DIR) -name '*.cpp') $(LDFLAGS) -o $@

# Run tests
test: $(TEST_EXES)
	@echo "Running tests..."
	@for test in $(TEST_EXES); do ./$$test || exit 1; done

# Clean up tests
clean:
	rm -rf $(TEST_EXES)

# Show configuration
.PHONY: info
info:
	@echo "Platform: $(PLATFORM)"
	@echo "Compiler: $(CXX)"
	@echo "CXXFLAGS: $(CXXFLAGS)"
	@echo "LDFLAGS: $(LDFLAGS)"
