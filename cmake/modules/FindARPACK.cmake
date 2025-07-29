# FindARPACK.cmake
# Find the ARPACK library
#
# The following variables are set if ARPACK is found:
#  ARPACK_FOUND        - True when ARPACK is found
#  ARPACK_INCLUDE_DIRS - The directory where ARPACK headers are located
#  ARPACK_LIBRARIES    - The libraries to link against for ARPACK
#
# The user can set the following variables to guide the search:
#  ARPACK_ROOT         - The root directory where ARPACK is installed
#  ARPACK_INCLUDEDIR   - The directory where ARPACK headers are installed
#  ARPACK_LIBRARYDIR   - The directory where ARPACK libraries are installed

# Look for the arpack header (Note: ARPACK doesn't typically have its own headers)
find_path(ARPACK_INCLUDE_DIR
  NAMES arpack.h arpack.hpp arpackdef.h
  HINTS
    ${ARPACK_INCLUDEDIR}
    ${ARPACK_ROOT}/include
    $ENV{ARPACK_ROOT}/include
  PATH_SUFFIXES
    arpack
)

# ARPACK might not have headers, so we don't require them
if(NOT ARPACK_INCLUDE_DIR)
  set(ARPACK_INCLUDE_DIR "")
  message(STATUS "ARPACK headers not found, but that's often normal")
endif()

# Look for the arpack library
find_library(ARPACK_LIBRARY
  NAMES arpack arpack-ng
  HINTS
    ${ARPACK_LIBRARYDIR}
    ${ARPACK_ROOT}/lib
    $ENV{ARPACK_ROOT}/lib
  PATH_SUFFIXES
    lib
    lib64
)

# Include dependent libraries
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

# Handle the QUIETLY and REQUIRED arguments
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(ARPACK
  DEFAULT_MSG
  ARPACK_LIBRARY
  BLAS_FOUND
  LAPACK_FOUND
)

# Set output variables
if(ARPACK_FOUND)
  set(ARPACK_LIBRARIES ${ARPACK_LIBRARY} ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
  set(ARPACK_INCLUDE_DIRS ${ARPACK_INCLUDE_DIR})
endif()

# Advanced variables for GUI users
mark_as_advanced(
  ARPACK_INCLUDE_DIR
  ARPACK_LIBRARY
) 