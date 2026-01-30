cmake_minimum_required(VERSION 3.15)
project(VegasAfterglow LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  add_compile_options(-O3 -DNDEBUG -pipe -w -DXTENSOR_DISABLE_EXCEPTIONS -march=native)
  add_compile_options(-fno-math-errno -fno-trapping-math -fno-signed-zeros -ffp-contract=fast)
elseif(MSVC)
  add_compile_options(/O2 /DNDEBUG /MP /GL /std:c++20 /DXTENSOR_DISABLE_EXCEPTIONS)
endif()

include(CheckIPOSupported)
check_ipo_supported(RESULT result OUTPUT output)
if(result)
  set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
  message(STATUS "LTO/IPO enabled")
else()
  message(WARNING "LTO/IPO is not supported: ${output}")
endif()

option(AFTERGLOW_PROFILE "Enable per-stage profiling instrumentation" OFF)

add_compile_definitions(
  NO_XTENSOR_IO
  $<$<BOOL:${AFTERGLOW_PROFILE}>:AFTERGLOW_PROFILE>
)

find_package(Python COMPONENTS Interpreter Development.Module NumPy REQUIRED)
find_package(pybind11 CONFIG REQUIRED)

if(APPLE)
    list(PREPEND CMAKE_PREFIX_PATH "/opt/homebrew/opt/libomp" "/usr/local/opt/libomp" "/opt/local")
endif()

find_package(OpenMP REQUIRED)
message(STATUS "OpenMP found - enabling parallel batch processing")

file(GLOB_RECURSE ALL_SOURCES
  "${CMAKE_CURRENT_SOURCE_DIR}/pybind/*.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp"
)

pybind11_add_module(VegasAfterglowC ${ALL_SOURCES})

target_link_libraries(VegasAfterglowC PUBLIC pybind11::module Python::NumPy)
target_link_libraries(VegasAfterglowC PUBLIC OpenMP::OpenMP_CXX)
target_compile_definitions(VegasAfterglowC PRIVATE HAVE_OPENMP)


target_include_directories(VegasAfterglowC PRIVATE
  "${CMAKE_CURRENT_SOURCE_DIR}/include"
  "${CMAKE_CURRENT_SOURCE_DIR}/src"
  "${CMAKE_CURRENT_SOURCE_DIR}/external"
)


set_target_properties(VegasAfterglowC PROPERTIES
  LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/VegasAfterglow"
)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()

# Install the module into the Python package directory
install(TARGETS VegasAfterglowC DESTINATION VegasAfterglow)


# ============================================================================
# Test/Benchmark Executables
# ============================================================================

option(BUILD_TESTS "Build test executables" ON)

if(BUILD_TESTS)
  # Common source files needed for executables
  set(VEGAS_CORE_SOURCES
          src/radiation/synchrotron.cpp
          src/radiation/power-law-syn.cpp
          src/radiation/inverse-compton.cpp
          src/radiation/smooth-power-law-syn.cpp
          src/radiation/prompt.cpp
          src/dynamics/shock.cpp
          src/core/mesh.cpp
          src/core/observer.cpp
          src/core/physics.cpp
  )

  # Demo executable
  add_executable(demo
          tests/demo/demo.cpp
          ${VEGAS_CORE_SOURCES}
  )

  target_include_directories(demo PRIVATE
          ${CMAKE_CURRENT_SOURCE_DIR}/include
          ${CMAKE_CURRENT_SOURCE_DIR}/src
          ${CMAKE_CURRENT_SOURCE_DIR}/external
  )

  target_compile_features(demo PRIVATE cxx_std_20)

  # Set output directory for test executables
  set_target_properties(demo PROPERTIES
          RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tests"
  )

  # Add profiling support for demo
  target_compile_options(demo PRIVATE
          $<$<CONFIG:Profiling>:-pg -g -fno-omit-frame-pointer>
          $<$<CONFIG:RelWithDebInfo>:-g -fno-omit-frame-pointer>
  )

  target_link_options(demo PRIVATE
          $<$<CONFIG:Profiling>:-pg>
  )

  # Benchmark executables (optional - only build if files exist)
  if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/tests/cpp-benchmark/benchmark.cpp")
    add_executable(benchmark
            tests/cpp-benchmark/benchmark.cpp
            ${VEGAS_CORE_SOURCES}
    )

    target_include_directories(benchmark PRIVATE
            ${CMAKE_CURRENT_SOURCE_DIR}/include
            ${CMAKE_CURRENT_SOURCE_DIR}/src
            ${CMAKE_CURRENT_SOURCE_DIR}/external
    )

    target_compile_features(benchmark PRIVATE cxx_std_20)

    set_target_properties(benchmark PROPERTIES
            RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tests"
    )

    target_compile_options(benchmark PRIVATE
            $<$<CONFIG:Profiling>:-pg -g -fno-omit-frame-pointer>
            $<$<CONFIG:RelWithDebInfo>:-g -fno-omit-frame-pointer>
    )

    target_link_options(benchmark PRIVATE
            $<$<CONFIG:Profiling>:-pg>
    )
  endif()

  # Profile benchmark executable (optional)
  if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/tests/cpp-benchmark/profile_benchmark.cpp")
    add_executable(profile_benchmark
            tests/cpp-benchmark/profile_benchmark.cpp
            ${VEGAS_CORE_SOURCES}
    )

    target_include_directories(profile_benchmark PRIVATE
            ${CMAKE_CURRENT_SOURCE_DIR}/include
            ${CMAKE_CURRENT_SOURCE_DIR}/src
            ${CMAKE_CURRENT_SOURCE_DIR}/external
    )

    target_compile_features(profile_benchmark PRIVATE cxx_std_20)

    set_target_properties(profile_benchmark PROPERTIES
            RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tests"
    )
  endif()
endif()
