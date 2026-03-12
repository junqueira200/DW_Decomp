#----------------------------------------------------------------
# Generated CMake target import file for configuration "Debug".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "::DwDecomp" for configuration "Debug"
set_property(TARGET ::DwDecomp APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
set_target_properties(::DwDecomp PROPERTIES
  IMPORTED_LOCATION_DEBUG "${_IMPORT_PREFIX}/lib64/libDwDecompd.so"
  IMPORTED_SONAME_DEBUG "libDwDecompd.so"
  )

list(APPEND _cmake_import_check_targets ::DwDecomp )
list(APPEND _cmake_import_check_files_for_::DwDecomp "${_IMPORT_PREFIX}/lib64/libDwDecompd.so" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
