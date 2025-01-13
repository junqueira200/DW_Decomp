cmake_minimum_required(VERSION 3.31)

project("DwDecomp"
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

find_package(Eigen3 3.3 REQUIRED)# NO_MODULE)
include_directories(${EIGEN3_INCLUDE_DIR})
find_package(GUROBI REQUIRED)
include_directories(${GUROBI_INCLUDE_DIRS})
link_directories($ENV{GUROBI_HOME}/lib)
find_package(Boost)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -g -pedantic -Wall -Wextra -fPIC")# -fsanitize=address")# -Wall")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D __FILENAME__='\"$$(subst  ${CMAKE_SOURCE_DIR}/,,$$(abspath $$<))\"'")

add_definitions(-DVAR_VECTOR_SANITY_CHECK=1)


add_library(${PROJECT_NAME} SHARED)

target_sources(${PROJECT_NAME}
    PRIVATE
        src/DW_Decomp.cpp
        src/BranchAndPrice.cpp
        src/SearchStrategy.cpp
        src/PrimalHeuristic.cpp
        src/Branch.cpp
)

target_link_libraries(${PROJECT_NAME} Eigen3::Eigen)
target_link_libraries(${PROJECT_NAME} gurobi_c++)
target_link_libraries(${PROJECT_NAME} gurobi)
target_link_libraries(${PROJECT_NAME} Boost::boost)


                        
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
    OUTPUT "${CMAKE_BINARY_DIR}/include-directories.txt"
    CONTENT ${includeDirectories}
)

# without it public headers won't get installed
set(public_headers
    include/DW_Decomp.h
    include/Aux.h
    include/Sparse.h
    include/SparseOp.h
    include/safe_matrix.h
    include/safe_vector.h
    include/BranchAndPrice.h
    include/SearchStrategy.h
    include/PrimalHeuristic.h
    include/Branch.h
    include/Alarm.h
    include/Statistics.h
)

# where to find our CMake modules
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
#message(STATUS "CMake module path: ${CMAKE_MODULE_PATH}")
include(Installing)
