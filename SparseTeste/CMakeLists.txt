cmake_minimum_required(VERSION 3.26)

project("SparseTeste")

set(CMAKE_CXX_STANDARD 23)

add_compile_definitions(EIGEN_DONT_PARALLELIZE)
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
find_package(Eigen3 3.3 REQUIRED NO_MODULE)
find_package(GUROBI REQUIRED)
include_directories(${GUROBI_INCLUDE_DIRS})
include_directories(include)
link_directories($ENV{GUROBI_HOME}/lib)


add_definitions(-DVAR_VECTOR_SANITY_CHECK=1)
message("set vector_sanity_check=1")

#add_library(DW_Decomp SHARED)
#find_package(DW_Decomp CONFIG REQUIRED)
#add_directory(../DW_Decomp)
#include_directories(../DW_Decomp/include)

#add_library(Grafo SHARED)
#find_package(Grafo CONFIG REQUIRED)
#add_directory(../Grafo)
#include_directories(../Grafo/include)

add_executable(${PROJECT_NAME})
target_sources(${PROJECT_NAME}
    PRIVATE
    src/main.cpp
)

target_link_libraries(${PROJECT_NAME} PRIVATE DwDecomp Grafo Eigen3::Eigen gurobi_c++ gurobi)

