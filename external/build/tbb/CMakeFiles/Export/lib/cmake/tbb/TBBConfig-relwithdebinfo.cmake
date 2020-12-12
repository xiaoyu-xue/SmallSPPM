#----------------------------------------------------------------
# Generated CMake target import file for configuration "RelWithDebInfo".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "TBB::tbb" for configuration "RelWithDebInfo"
set_property(TARGET TBB::tbb APPEND PROPERTY IMPORTED_CONFIGURATIONS RELWITHDEBINFO)
set_target_properties(TBB::tbb PROPERTIES
  IMPORTED_IMPLIB_RELWITHDEBINFO "${_IMPORT_PREFIX}/lib/tbb.lib"
  IMPORTED_LOCATION_RELWITHDEBINFO "${_IMPORT_PREFIX}/bin/tbb.dll"
  )

list(APPEND _IMPORT_CHECK_TARGETS TBB::tbb )
list(APPEND _IMPORT_CHECK_FILES_FOR_TBB::tbb "${_IMPORT_PREFIX}/lib/tbb.lib" "${_IMPORT_PREFIX}/bin/tbb.dll" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
