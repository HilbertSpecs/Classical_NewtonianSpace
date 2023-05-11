cmake_minimum_required(VERSION 3.20.0)
project(FindNewtonianSpacePackage CXX)

find_package(NewtonianSpace REQUIRED)

include(CMakePrintHelpers)
message("CMAKE_PREFIX_PATH: ${CMAKE_PREFIX_PATH}")
message("NEWTONIANSPACE_FOUND: ${NewtonianSpace_FOUND}")
#[[cmake_print_properties(TARGETS "NewtonianSpace::newtonianspace_shared" PROPERTIES
  IMPORTED_CONFIGURATIONS
  INCLUDE_DIRECTORIES
  INTERFACE_INCLUDE_DIRECTORIES
)]]
cmake_print_properties(TARGETS "newtonianspace_shared" PROPERTIES
  IMPORTED_CONFIGURATIONS
  INCLUDE_DIRECTORIES
  INTERFACE_INCLUDE_DIRECTORIES
)

add_executable(test main_ElementaryAnalysis.cpp)
#target_link_libraries(test NewtonianSpace::newtonianspace_shared)
target_link_libraries(test newtonianspace_shared)