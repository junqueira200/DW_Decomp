cmake_minimum_required(VERSION 3.5)
project(DW_Decomp)

set(CMAKE_CXX_STANDARD 20)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR})

find_package(Eigen3 3.3 REQUIRED NO_MODULE)
find_package(GUROBI REQUIRED)
include_directories(${GUROBI_INCLUDE_DIRS})
link_directories($ENV{GUROBI_HOME}/lib)

set(CMAKE_CXX_FLAGS "-O0 -g3")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D __FILENAME__='\"$(subst  ${CMAKE_SOURCE_DIR}/,,$(abspath $<))\"'")

add_definitions(-DVAR_VECTOR_SANITY_CHECK=1)

add_executable(dw_decomp main.cpp
        MNFP_Inst.cpp
        MNFP_Inst.h
        DW_Decomp.cpp
        DW_Decomp.h
        Aux.h)
target_link_libraries(dw_decomp Eigen3::Eigen)
target_link_libraries(dw_decomp gurobi_c++)
target_link_libraries(dw_decomp gurobi)

