cmake_minimum_required(VERSION 3.5)

project("Grafo"
    VERSION 0
    DESCRIPTION "Some library"
)

if (NOT CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    message(STATUS "This project has a top-level one called [${CMAKE_PROJECT_NAME}]")
else()
    message(STATUS "This project is a top-level one")
endif()

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

message(CMAKE_MODULE_PATH="${CMAKE_SOURCE_DIR}")

#set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

find_package(Eigen3 3.3 REQUIRED)# NO_MODULE)
include_directories(${EIGEN3_INCLUDE_DIR})

#set(CMAKE_CXX_FLAGS "-O0 -g3")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pedantic -O3")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D __FILENAME__='\"$$(subst  ${CMAKE_SOURCE_DIR}/,,$$(abspath $$<))\"'")



add_library(${PROJECT_NAME} SHARED)

target_sources(${PROJECT_NAME}
    PRIVATE
        src/Grafo.cpp
)

target_link_libraries(${PROJECT_NAME} Eigen3::Eigen)

                        
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
    OUTPUT "${CMAKE_BINARY_DIR}/include-directories-Grafo.txt"
    CONTENT ${includeDirectories}
)

# without it public headers won't get installed
set(public_headers
    include/GAux.h
    include/Grafo.h
)

# not for MSVC
if(CMAKE_COMPILER_IS_GNUCXX)
    # compile options for this target only
    # treat warnings as errors
    target_compile_options(${PROJECT_NAME} PRIVATE -Werror)
endif()

# where to find our CMake modules
#set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
#message(STATUS "CMake module path: ${CMAKE_MODULE_PATH}")
#include(Installing)
