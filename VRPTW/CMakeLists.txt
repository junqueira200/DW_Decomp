cmake_minimum_required(VERSION 3.31)

project("Vrptw")

set(CMAKE_CXX_STANDARD 23)
#set(CMAKE_CXX_SCAN_FOR_MODULES ON)

set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
find_package(Eigen3 3.3 REQUIRED NO_MODULE)
find_package(GUROBI REQUIRED)
include_directories(${GUROBI_INCLUDE_DIRS})
include_directories(include)
link_directories($ENV{GUROBI_HOME}/lib)


add_definitions(-DVAR_VECTOR_SANITY_CHECK=1)
message("set vector_sanity_check=1")


#set(CMAKE_CXX_FLAGS "-pedantic")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -g -fPIC")# -pedantic -Wall -Wextra -fPIC")

find_package(Boost CONFIG)

#add_library(DW_Decomp SHARED)
#find_package(DW_Decomp CONFIG REQUIRED)
#add_directory(../DW_Decomp)
#include_directories(../DW_Decomp/include)

#add_library(Grafo SHARED)
#find_package(Grafo CONFIG REQUIRED)
#add_directory(../Grafo)
#include_directories(../Grafo/include)

add_executable(${PROJECT_NAME}
    include/LabelingConstants.h
    include/Label.h
    include/NgSet.h
    include/LabelingUtils.h
    src/Label.cpp
    src/LabelingUtils.cpp
    include/Test.h
    src/Test.cpp)


#target_sources(${PROJECT_NAME}
        #    PUBLIC
        #        FILE_SET cxx_modules TYPE CXX_MODULES
#        FILES
        #            src/teste.ixx

        #    PUBLIC
        #        src/teste.cxx
#)

target_sources(${PROJECT_NAME}
    PRIVATE
        src/LabelingAlgorithm.cpp
        src/VrpTW_DecompLabeling.cpp
        src/VrpTW_Decomp.cpp
        src/Instancia.cpp
        src/SolomonInst.cpp
        src/LowerBound.cpp
        src/Label.cpp
        src/main.cpp
)

target_link_libraries(${PROJECT_NAME} PRIVATE DwDecomp Eigen3::Eigen gurobi_c++ gurobi Boost::boost Grafo CvrpSep)

