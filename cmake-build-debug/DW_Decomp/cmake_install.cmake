# Install script for directory: /Users/igor/Documents/[08] Projetos/[16]_DW_Decomp/DW_Decomp

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/Users/igor/Documents/[08] Projetos/[16]_DW_Decomp/install")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set path to fallback-tool for dependency-resolution.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/usr/bin/objdump")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES "/Users/igor/Documents/[08] Projetos/[16]_DW_Decomp/cmake-build-debug/DW_Decomp/libDwDecomp.dylib")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libDwDecomp.dylib" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libDwDecomp.dylib")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Library/gurobi1203/macos_universal2/lib"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libDwDecomp.dylib")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" -x "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libDwDecomp.dylib")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/DwDecomp" TYPE FILE FILES
    "/Users/igor/Documents/[08] Projetos/[16]_DW_Decomp/DW_Decomp/include/DW_Decomp.h"
    "/Users/igor/Documents/[08] Projetos/[16]_DW_Decomp/DW_Decomp/include/Aux.h"
    "/Users/igor/Documents/[08] Projetos/[16]_DW_Decomp/DW_Decomp/include/Sparse.h"
    "/Users/igor/Documents/[08] Projetos/[16]_DW_Decomp/DW_Decomp/include/SparseOp.h"
    "/Users/igor/Documents/[08] Projetos/[16]_DW_Decomp/DW_Decomp/include/safe_matrix.h"
    "/Users/igor/Documents/[08] Projetos/[16]_DW_Decomp/DW_Decomp/include/safe_vector.h"
    "/Users/igor/Documents/[08] Projetos/[16]_DW_Decomp/DW_Decomp/include/safe_3D_matrix.h"
    "/Users/igor/Documents/[08] Projetos/[16]_DW_Decomp/DW_Decomp/include/BranchAndPrice.h"
    "/Users/igor/Documents/[08] Projetos/[16]_DW_Decomp/DW_Decomp/include/SearchStrategy.h"
    "/Users/igor/Documents/[08] Projetos/[16]_DW_Decomp/DW_Decomp/include/PrimalHeuristic.h"
    "/Users/igor/Documents/[08] Projetos/[16]_DW_Decomp/DW_Decomp/include/Branch.h"
    "/Users/igor/Documents/[08] Projetos/[16]_DW_Decomp/DW_Decomp/include/Alarm.h"
    "/Users/igor/Documents/[08] Projetos/[16]_DW_Decomp/DW_Decomp/include/Statistics.h"
    )
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/cmake/DwDecompTargets.cmake")
    file(DIFFERENT _cmake_export_file_changed FILES
         "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/cmake/DwDecompTargets.cmake"
         "/Users/igor/Documents/[08] Projetos/[16]_DW_Decomp/cmake-build-debug/DW_Decomp/CMakeFiles/Export/272ceadb8458515b2ae4b5630a6029cc/DwDecompTargets.cmake")
    if(_cmake_export_file_changed)
      file(GLOB _cmake_old_config_files "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/cmake/DwDecompTargets-*.cmake")
      if(_cmake_old_config_files)
        string(REPLACE ";" ", " _cmake_old_config_files_text "${_cmake_old_config_files}")
        message(STATUS "Old export file \"$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/cmake/DwDecompTargets.cmake\" will be replaced.  Removing files [${_cmake_old_config_files_text}].")
        unset(_cmake_old_config_files_text)
        file(REMOVE ${_cmake_old_config_files})
      endif()
      unset(_cmake_old_config_files)
    endif()
    unset(_cmake_export_file_changed)
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/cmake" TYPE FILE FILES "/Users/igor/Documents/[08] Projetos/[16]_DW_Decomp/cmake-build-debug/DW_Decomp/CMakeFiles/Export/272ceadb8458515b2ae4b5630a6029cc/DwDecompTargets.cmake")
  if(CMAKE_INSTALL_CONFIG_NAME MATCHES "^()$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/cmake" TYPE FILE FILES "/Users/igor/Documents/[08] Projetos/[16]_DW_Decomp/cmake-build-debug/DW_Decomp/CMakeFiles/Export/272ceadb8458515b2ae4b5630a6029cc/DwDecompTargets-noconfig.cmake")
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/cmake" TYPE FILE FILES
    "/Users/igor/Documents/[08] Projetos/[16]_DW_Decomp/cmake-build-debug/DW_Decomp/DwDecompConfig.cmake"
    "/Users/igor/Documents/[08] Projetos/[16]_DW_Decomp/cmake-build-debug/DW_Decomp/DwDecompConfigVersion.cmake"
    )
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
if(CMAKE_INSTALL_LOCAL_ONLY)
  file(WRITE "/Users/igor/Documents/[08] Projetos/[16]_DW_Decomp/cmake-build-debug/DW_Decomp/install_local_manifest.txt"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
endif()
