cmake_minimum_required(VERSION 3.15)
project(VegasAfterglow LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
option(AFTERGLOW_ENABLE_WARNINGS "Enable compiler warnings" OFF)

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  add_compile_options(-O3 -pipe -DNDEBUG -DXTENSOR_DISABLE_EXCEPTIONS -DXTENSOR_USE_XSIMD)
  add_compile_options(-fno-math-errno -fno-trapping-math -fno-signed-zeros -ffp-contract=fast -w)
elseif(MSVC)
  add_compile_options(/O2 /DNDEBUG /MP /GL /std:c++20 /DXTENSOR_DISABLE_EXCEPTIONS /DXTENSOR_USE_XSIMD)
endif()

if(AFTERGLOW_ENABLE_WARNINGS)
  if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    add_compile_options(-Wall -Wextra -Wpedantic)
  elseif(MSVC)
    add_compile_options(/W4)
  endif()
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

if(APPLE AND CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    # AppleClang doesn't have built-in OpenMP; use Homebrew's libomp (static)
    foreach(_omp_root IN ITEMS "/opt/homebrew/opt/libomp" "/usr/local/opt/libomp")
        if(EXISTS "${_omp_root}/lib/libomp.a")
            set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp" CACHE STRING "" FORCE)
            set(OpenMP_CXX_LIB_NAMES "omp" CACHE STRING "" FORCE)
            set(OpenMP_omp_LIBRARY "${_omp_root}/lib/libomp.a" CACHE FILEPATH "" FORCE)
            include_directories("${_omp_root}/include")
            message(STATUS "Using static OpenMP: ${_omp_root}/lib/libomp.a")
            break()
        endif()
    endforeach()
endif()

find_package(OpenMP REQUIRED)
message(STATUS "OpenMP found - enabling parallel batch processing")

set(AFTERGLOW_SOURCES
  pybind/mcmc.cpp
  pybind/pybind.cpp
  pybind/pymodel.cpp
  src/core/mesh.cpp
  src/core/observer.cpp
  src/core/physics.cpp
  src/dynamics/shock.cpp
  src/radiation/inverse-compton.cpp
  src/radiation/power-law-syn.cpp
  src/radiation/prompt.cpp
  src/radiation/smooth-power-law-syn.cpp
  src/radiation/synchrotron.cpp
  src/util/IO.cpp
  src/util/utilities.cpp
)

pybind11_add_module(VegasAfterglowC ${AFTERGLOW_SOURCES})

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
