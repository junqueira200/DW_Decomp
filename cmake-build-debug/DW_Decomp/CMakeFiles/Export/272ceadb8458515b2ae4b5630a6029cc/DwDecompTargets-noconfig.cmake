#----------------------------------------------------------------
# Generated CMake target import file.
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "::DwDecomp" for configuration ""
set_property(TARGET ::DwDecomp APPEND PROPERTY IMPORTED_CONFIGURATIONS NOCONFIG)
set_target_properties(::DwDecomp PROPERTIES
  IMPORTED_LOCATION_NOCONFIG "${_IMPORT_PREFIX}/lib/libDwDecomp.dylib"
  IMPORTED_SONAME_NOCONFIG "@rpath/libDwDecomp.dylib"
  )

list(APPEND _cmake_import_check_targets ::DwDecomp )
list(APPEND _cmake_import_check_files_for_::DwDecomp "${_IMPORT_PREFIX}/lib/libDwDecomp.dylib" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
