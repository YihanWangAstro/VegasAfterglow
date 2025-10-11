cmake_minimum_required(VERSION 3.15)
project(VegasAfterglow LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  add_compile_options(-O3 -flto -DNDEBUG -pipe -w -DXTENSOR_DISABLE_EXCEPTIONS -DXTENSOR_USE_XSIMD)
elseif(MSVC)
  add_compile_options(/O2 /DNDEBUG /MP /GL /std:c++20 /DXTENSOR_DISABLE_EXCEPTIONS /DXTENSOR_USE_XSIMD)
  add_link_options(/LTCG)
endif()

add_compile_definitions(
  NO_XTENSOR_IO
)

find_package(Python COMPONENTS Interpreter Development.Module NumPy REQUIRED)
find_package(pybind11 CONFIG REQUIRED)

file(GLOB_RECURSE ALL_SOURCES
  "${CMAKE_CURRENT_SOURCE_DIR}/pybind/*.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp"
)

pybind11_add_module(VegasAfterglowC ${ALL_SOURCES})

target_link_libraries(VegasAfterglowC PUBLIC pybind11::module Python::NumPy)

target_include_directories(VegasAfterglowC PRIVATE
  "${CMAKE_CURRENT_SOURCE_DIR}/include"
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
