cmake_minimum_required(VERSION 3.5)
project(BA)
#set(gtest_build_samples 1)
#add_subdirectory(googletest)

#set(BUILD_SHARED_LIBS OFF)
#SET(BUILD_SHARED_LIBRARIES OFF)
#set(ITK_DYNAMIC_LOADING OFF)
#set(CMAKE_EXE_LINKER_FLAGS "-static")
#set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")


#set(CPLEX_FIND_QUIETLY OFF)
#set(CGAL_FIND_QUIETLY OFF)
set(UTIL_BOOST_LIBRARIES system filesystem program_options thread )
if(NOT ALG_CMAKE_DIR)
    if(NOT "$ENV{ALG_CMAKE_DIR}" STREQUAL "")
        set(ALG_CMAKE_DIR "$ENV{ALG_CMAKE_DIR}")
    else()
        set(ALG_CMAKE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
    endif()
endif()


include(${ALG_CMAKE_DIR}/util.cmake)
#set(CPLEX_ROOT /opt/ibm/ILOG/CPLEX_Studio1271/)


find_package(CGAL REQUIRED)
set(Boost_USE_STATIC_LIBS   ON)
find_package(Boost COMPONENTS system filesystem program_options thread REQUIRED)

find_package(CPLEX REQUIRED)



set(ALGCPLEX_INCLUDE_DIR ${UTIL_ALGCPLEX_ROOT}/include/algcplex)
set(LIBRARIES ${LIBRARIES} ${CPLEX_LIBRARIES} ${CGAL_LIBRARIES} ${Boost_LIBRARIES} pthread)
set(INCLUDE_DIRS ${INCLUDE_DIRS} ${CPLEX_INCLUDE_DIRS} ${ALGCPLEX_INCLUDE_DIR} ${CGAL_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS})

set(SOURCE_FILES Capturing/Capturing2D.cpp Capturing/Header/Capturing2D.h Capturing/Segment.cpp Capturing/Header/Segment.h Utilities/GeoUtil.cpp Utilities/GeoUtil.h Utilities/Graph/undirectedMatrixGraph.cpp Utilities/Graph/undirectedMatrixGraph.h Utilities/Graph/undirectedAdjGraph.cpp Utilities/Graph/undirectedAdjGraph.h Utilities/Graph/undirectedGraph.cpp Utilities/Graph/undirectedGraph.h Visualization/SegmentVisualization.cpp Visualization/SegmentVisualization.h Visualization/simple_svg_1.0.0.hpp Capturing/CapturingGraph.cpp Capturing/Header/CapturingGraph.h Capturing/CapturingSolution.cpp Capturing/Header/CapturingSolution.h Visualization/Imports.cpp Visualization/Imports.h Capturing/Generator/Cap2DGenerator.cpp Capturing/Generator/Header/Cap2DGenerator.h MacroTestingFramework/MacroTestingInstance.hpp MacroTestingFramework/TestInstance.hpp "MacroTestingFramework/Wrapper Templates/Generator.h" "MacroTestingFramework/Wrapper Templates/Solver.h" "MacroTestingFramework/Wrapper Templates/Solution.h" Capturing/Wrapper/GeneratorWrapper.hpp Capturing/Wrapper/SolverWrapper.hpp MacroTestingFramework/TestInstanceContainer.h Capturing/Metaheuristics/SimulatedAnnealing.cpp Capturing/Metaheuristics/SimulatedAnnealing.h Capturing/Metaheuristics/EvolutionComputation.cpp Capturing/Metaheuristics/EvolutionComputation.h)
add_executable(BA main.cpp ${SOURCE_FILES})
util_setup_target(BA)

# add other libraries (CPLEX)
target_include_directories(BA PRIVATE ${INCLUDE_DIRS})
target_link_libraries(BA ${LIBRARIES})
message(${LIBRARIES})
#add_executable(Visualization Visualization/main_1.0.0.cpp Visualization/simple_svg_1.0.0.hpp)
#util_setup_target(Visualization)

# add other libraries (CPLEX)
#target_include_directories(Visualization PRIVATE ${INCLUDE_DIRS})
#target_link_libraries(Visualization ${LIBRARIES})

set(TEST_SOURCE_FILES Tests/GraphTests.cpp Tests/MacroTestFrameWorkTests.cpp Capturing/Metaheuristics/TabuSearch.cpp Capturing/Metaheuristics/TabuSearch.h Tests/CGALTests.cpp)
add_executable(testGoogle ${TEST_SOURCE_FILES} ${SOURCE_FILES})
#include_directories(./googletest )
#link_directories(./googletest)
target_include_directories(testGoogle PRIVATE gtest gtest_main ${INCLUDE_DIRS})
target_link_libraries(testGoogle gtest gtest_main ${LIBRARIES})
