cmake_minimum_required(VERSION 3.26)

project("2L_SDVRP"
    VERSION 0
    DESCRIPTION "Some library"
)

if (NOT CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    message(STATUS "This project has a top-level one called [${CMAKE_PROJECT_NAME}]")
else()
    message(STATUS "This project is a top-level one")
endif()

set(CMAKE_CXX_STANDARD 23)
add_definitions(-DVAR_VECTOR_SANITY_CHECK=1)
message("set vector_sanity_check=1")


message(CMAKE_MODULE_PATH="${CMAKE_SOURCE_DIR}")

set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

#find_package(Eigen3 3.3 REQUIRED)# NO_MODULE)
include_directories(${EIGEN3_INCLUDE_DIR})
find_package(ortools CONFIG REQUIRED)
#find_package(GUROBI REQUIRED)
#include_directories(${GUROBI_INCLUDE_DIRS})
#link_directories($ENV{GUROBI_HOME}/lib)
#find_package(Boost)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -g -pedantic -fPIC")# -Wall -Wextra -fPIC")# -fsanitize=address")# -Wall")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D __FILENAME__='\"$$(subst  ${CMAKE_SOURCE_DIR}/,,$$(abspath $$<))\"'")

add_definitions(-DVAR_VECTOR_SANITY_CHECK=1)
add_definitions(-DVAR_SOLUTION_BORROW_CHECKER=1)
message("set solution_borrow_checker=1")

add_library(${PROJECT_NAME} SHARED)

target_sources(${PROJECT_NAME}
    PRIVATE
        src/BuscaLocal.cpp
        src/Construtivo.cpp
        src/rand.cpp
        src/InputOutput.cpp
        src/Instancia.cpp
        src/Solucao.cpp
        src/ConstrutivoBin.cpp
        src/Ig.cpp
        src/BinPackingCP.cpp
)

target_link_libraries(${PROJECT_NAME} Eigen3::Eigen)
target_link_libraries(${PROJECT_NAME} ortools::ortools)
#target_link_libraries(${PROJECT_NAME} DwDecomp)
#target_link_libraries(${PROJECT_NAME} gurobi_c++)
#target_link_libraries(${PROJECT_NAME} gurobi)
#target_link_libraries(${PROJECT_NAME} Boost::boost)


                        
target_include_directories(${PROJECT_NAME}
    PRIVATE
        # where the library itself will look for its internal headers
        ${CMAKE_CURRENT_SOURCE_DIR}/src
    PUBLIC
        # where top-level project will look for the library's public headers
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        # where external projects will look for the library's public headers
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

# totally optional listing of include directories
get_property(inclds
    TARGET ${PROJECT_NAME}
    PROPERTY INCLUDE_DIRECTORIES
)
message(STATUS "Listing include directories for ${PROJECT_NAME}...")
# "dumb" listing with "raw" generator expressions on configuration
foreach(pth ${inclds})
    message(STATUS "  ${pth}")
endforeach()
# actually evaluated generator expressions printed to file on generation
string(JOIN "\n" includeDirectories ${inclds})
file(GENERATE
    OUTPUT "${CMAKE_BINARY_DIR}/include-directories-2L_SDVRP.txt"
    CONTENT ${includeDirectories}
)

# without it public headers won't get installed
set(public_headers
    include/AuxT.h
    include/BinPackingCP.h
    include/BuscaLocal.h
    include/Constantes.h
    include/ConstrutivoBin.h
    include/Construtivo.h
    include/Convencoes.h
    include/cxxopts.hpp
    include/Ig.h
    include/InputOutput.h
    include/Instancia.h
    include/rand.h
    include/safe_matrix.h
    include/safe_vector.h
    include/sefe_array.h
    include/Solucao.h
)

# where to find our CMake modules
#set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
#message(STATUS "CMake module path: ${CMAKE_MODULE_PATH}")
#include(Installing)
