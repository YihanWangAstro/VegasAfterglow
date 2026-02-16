cmake_minimum_required(VERSION 3.15)
project(VegasAfterglow LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  add_compile_options(-O3 -pipe -DNDEBUG -DXTENSOR_DISABLE_EXCEPTIONS -Wall -Wextra -Wpedantic)
  add_compile_options(-fno-math-errno -fno-trapping-math -fno-signed-zeros -ffp-contract=fast
                       -freciprocal-math)
  add_compile_options(-fvisibility=hidden)
  if(APPLE)
    add_compile_options(-fveclib=Accelerate)
  endif()
elseif(MSVC)
  add_compile_options(/O2 /DNDEBUG /MP /GL /std:c++20 /DXTENSOR_DISABLE_EXCEPTIONS)
endif()

include(CheckIPOSupported)
check_ipo_supported(RESULT result OUTPUT output)
if(result)
  set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
  if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    add_compile_options(-flto=full)
    add_link_options(-flto=full)
  endif()
  message(STATUS "LTO/IPO enabled (full)")
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

set(AFTERGLOW_SOURCES
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


target_include_directories(VegasAfterglowC PRIVATE
  "${CMAKE_CURRENT_SOURCE_DIR}/include"
  "${CMAKE_CURRENT_SOURCE_DIR}/src"
)
target_include_directories(VegasAfterglowC SYSTEM PRIVATE
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
