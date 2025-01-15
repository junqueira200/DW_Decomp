cmake_minimum_required(VERSION 3.31)

project("CvrpSep"
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


set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -g  -fPIC")# -fsanitize=address")# -Wall")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D __FILENAME__='\"$$(subst  ${CMAKE_SOURCE_DIR}/,,$$(abspath $$<))\"'")

add_definitions(-DVAR_VECTOR_SANITY_CHECK=1)


add_library(${PROJECT_NAME} SHARED)

target_sources(${PROJECT_NAME}
    PRIVATE
   src/basegrph.c
   src/binpack.c
   src/blocks.c
   src/brnching.c
   src/capsep.c
   src/cnstrmgr.c
   src/combsep.c
   src/compcuts.c
   src/compress.c
   src/cutbase.c
   src/fcapfix.c
   src/fcisep.c
   src/fcits.c
   src/glmsep.c
   src/grsearch.c
   src/hpmstar.c
   src/htoursep.c
   src/intap.c
   src/memmod.c
   src/mstarsep.c
   src/mxf.c
   src/newhtour.c
   src/sort.c
   src/strcomb.c
   src/strngcmp.c
   src/twomatch.c
)

                        
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
    OUTPUT "${CMAKE_BINARY_DIR}/include-directories-CvrpSep.txt"
    CONTENT ${includeDirectories}
)

# without it public headers won't get installed
set(public_headers
    include/basegrph.h
    include/binpack.h
    include/blocks.h
    include/brnching.h
    include/capsep.h
    include/cnstrmgr.h
    include/combsep.h
    include/compcuts.h
    include/compress.h
    include/cutbase.h
    include/fcapfix.h
    include/fcisep.h
    include/fcits.h
    include/glmsep.h
    include/grsearch.h
    include/hpmstar.h
    include/htoursep.h
    include/intap.h
    include/memmod.h
    include/mstarsep.h
    include/mxf.h
    include/newhtour.h
    include/sort.h
    include/strcomb.h
    include/strngcmp.h
    include/twomatch.h
)

# where to find our CMake modules
#set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
#message(STATUS "CMake module path: ${CMAKE_MODULE_PATH}")
#include(Installing)
