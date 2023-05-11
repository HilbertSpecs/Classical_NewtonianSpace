add_library(newtonianspace_obj OBJECT NewtonianSpace.cpp)
target_include_directories(newtonianspace_obj INTERFACE
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
  "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
)
set_target_properties(newtonianspace_obj PROPERTIES
    PUBLIC_HEADER src/newtonianspace/include/newtonianspace/NewtonianSpace.hpp
    POSITION_INDEPENDENT_CODE 1
)

add_library(newtonianspace_shared SHARED)
target_link_libraries(newtonianspace_shared newtonianspace_obj)
add_library(newtonianspace_static STATIC)
target_link_libraries(newtonianspace_static newtonianspace_obj)

#Add this to create executable target in full_project
#[[add_executable(elementaryanalysis ElementaryAnalysis.cpp)
target_link_libraries(elementaryanalysis newtonianspace_static)]]
#EndAdd

include(Coverage)
CleanCoverage(newtonianspace_static)
include(Format)
Format(newtonianspace_static .)
include(CppCheck)
AddCppCheck(newtonianspace_obj)
include(Doxygen)
Doxygen(newtonianspace src/newtonianspace)