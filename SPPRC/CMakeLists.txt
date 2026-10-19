cmake_minimum_required(VERSION 3.16)

project("SPPRC"
    VERSION 0
    DESCRIPTION "Some library"
)

if (NOT CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    message(STATUS "This project has a top-level one called [${CMAKE_PROJECT_NAME}]")
else()
    message(STATUS "This project is a top-level one")
endif()

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_BUILD_TYPE Release)
add_definitions(-DVAR_VECTOR_SANITY_CHECK=1)
message("set vector_sanity_check=1")

message(CMAKE_MODULE_PATH="${CMAKE_SOURCE_DIR}")

set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

cmake_policy(SET CMP0074 NEW)
cmake_policy(SET CMP0167 NEW)
cmake_policy(SET CMP0144 NEW)
find_package(Eigen3 3.3 REQUIRED)# NO_MODULE)
include_directories(${EIGEN3_INCLUDE_DIR})
#find_package(ortools CONFIG REQUIRED)
#find_package(GUROBI REQUIRED)
#include_directories(${GUROBI_INCLUDE_DIRS})
#link_directories($ENV{GUROBI_HOME}/lib)

#include_directories($ENV{SCIPpp}/include)
#link_directories($ENV{SCIPpp}/lib)


find_package(Boost)
find_package(OpenMP REQUIRED)


set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -g -fPIC")# -Wall -Wextra -fPIC")# -fsanitize=address")# -Wall")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D __FILENAME__='\"$$(subst  ${CMAKE_SOURCE_DIR}/,,$$(abspath $$<))\"'")

#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -ftrivial-auto-var-init=zero")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ftrivial-auto-var-init=zero")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address,undefined")

add_definitions(-DVAR_VECTOR_SANITY_CHECK=1)
add_definitions(-DVAR_SOLUTION_BORROW_CHECKER=1)
message("set solution_borrow_checker=1")

add_library(${PROJECT_NAME} SHARED)

target_sources(${PROJECT_NAME}
    PRIVATE
    src/InstanceVRPTW.cpp
    src/Cvrp_DecompLabeling.cpp
    src/LabelingUtils.cpp
    src/VrpTW_DecompLabeling.cpp
    src/LowerBound.cpp
    src/Label.cpp
    src/SolomonInst.cpp
    src/LabelingAlgorithm.cpp
    src/VrpTW_Decomp.cpp
	src/Grafo.cpp
    src/rand.cpp
)


target_link_libraries(${PROJECT_NAME} PRIVATE Eigen3::Eigen)
target_link_libraries(${PROJECT_NAME} PRIVATE Boost::boost)
target_link_libraries(${PROJECT_NAME} PRIVATE OpenMP::OpenMP_CXX)

                        
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
    OUTPUT "${CMAKE_BINARY_DIR}/include-directories-SPPRC.txt"
    CONTENT ${includeDirectories}
)

# without it public headers won't get installed
set(public_headers
    include/Cvrp_DecompLabeling.h
    include/LabelingConstants.h
    include/MemoryPool.h
    include/Test.h
    include/Label.h
    include/LabelingUtils.h
    include/NgSet.h
    include/VrpTW_Decomp.h
    include/LabelingAlgorithm.h
    include/LowerBound.h
    include/SolomonInst.h
    include/VrpTW_DecompLabeling.h
	include/Aux.h
	include/safe_matrix.h
	include/safe_vector.h
	include/safe_3D_matrix.h
	include/Grafo.h
    include/InstanceVRPTW.h
    include/sefe_array.h
    include/rand.h
)

# where to find our CMake modules
#set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
#message(STATUS "CMake module path: ${CMAKE_MODULE_PATH}")
#include(Installing)
