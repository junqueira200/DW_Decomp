cmake_minimum_required(VERSION 3.31)

project("DW")

set(CMAKE_CXX_STANDARD 23)


set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
#find_package(Eigen3 3.3 REQUIRED NO_MODULE)
#find_package(GUROBI REQUIRED)
#include_directories(${GUROBI_INCLUDE_DIRS})
#link_directories($ENV{GUROBI_HOME}/lib)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -g -pedantic -Wall -Wextra -fPIC -march=native")#-Wall")

add_subdirectory(DW_Decomp)
include_directories(DW_Decomp/include)

add_subdirectory(CvrpSep)
include_directories(DW_Decomp/include)

add_subdirectory(Grafo)
include_directories(Grafo/include)

#add_subdirectory(Mnfp)
#include_directories(Mnfp/include)

add_subdirectory(SparseTeste)
include_directories(SparseTeste/include)

add_subdirectory(VRPTW)
include_directories(VRPTW/include)

#add_executable(${PROJECT_NAME})
#target_sources(${PROJECT_NAME}
#    PRIVATE
#    MNFP_Inst.cpp
#    main.cpp
#)

#target_link_libraries(${PROJECT_NAME} PRIVATE DwDecomp Grafo Eigen3::Eigen gurobi_c++ gurobi)

