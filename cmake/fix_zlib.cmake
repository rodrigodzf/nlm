# This file completely replaces the zlib CMakeLists.txt to prevent the shared library from being built
# and also patches matio's thirdParties.cmake to use only the static zlib
# Note: Static library configuration is handled globally by StaticLibraryConfig.cmake

# Wait until zlib is downloaded before patching
if(NOT EXISTS "${zlib_SOURCE_DIR}")
  return()
endif()

# Step 1: Patch zlib's CMakeLists.txt to only build the static library
# Read the original zlib CMakeLists.txt
file(READ "${zlib_SOURCE_DIR}/CMakeLists.txt" ZLIB_ORIGINAL_CMAKE)

# Find the line where the library definitions start
string(FIND "${ZLIB_ORIGINAL_CMAKE}" "#============================================================================\n# zlib" ZLIB_START_POS)
if(ZLIB_START_POS EQUAL -1)
  message(FATAL_ERROR "Could not find zlib library definition section in CMakeLists.txt")
endif()

# Extract the part of the file before the library definitions
string(SUBSTRING "${ZLIB_ORIGINAL_CMAKE}" 0 ${ZLIB_START_POS} ZLIB_HEADER)

# Configure our custom CMakeLists.txt with only the static library
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/cmake/patch/zlib_CMakeLists.txt.in" 
              "${CMAKE_CURRENT_BINARY_DIR}/zlib_CMakeLists.txt.patched" @ONLY)

# Read our patched file 
file(READ "${CMAKE_CURRENT_BINARY_DIR}/zlib_CMakeLists.txt.patched" ZLIB_PATCHED_CMAKE)

# Combine the original header with our patched library section
file(WRITE "${zlib_SOURCE_DIR}/CMakeLists.txt" "${ZLIB_HEADER}${ZLIB_PATCHED_CMAKE}")

message(STATUS "Successfully patched zlib to build only static libraries")

# Step 2: Patch matio's thirdParties.cmake to use zlibstatic
# Wait for matio to be downloaded
if(NOT EXISTS "${matio_SOURCE_DIR}")
  return()
endif()

# Set static libraries in matio's CMakeLists.txt directly
if(EXISTS "${matio_SOURCE_DIR}/CMakeLists.txt")
  file(READ "${matio_SOURCE_DIR}/CMakeLists.txt" MATIO_CMAKE)
  # Note: Static library configuration is now handled globally by StaticLibraryConfig.cmake
  # string(REPLACE "cmake_minimum_required" "set(BUILD_SHARED_LIBS OFF CACHE BOOL \"Build shared libraries\" FORCE)\ncmake_minimum_required" MATIO_CMAKE "${MATIO_CMAKE}")
  file(WRITE "${matio_SOURCE_DIR}/CMakeLists.txt" "${MATIO_CMAKE}")
  message(STATUS "Successfully patched matio's CMakeLists.txt to use only static libraries")
endif()

# Now patch the matio thirdParties.cmake
if(EXISTS "${matio_SOURCE_DIR}/cmake/thirdParties.cmake")
  # Find the MATIO_CREATE_ZLIB macro in thirdParties.cmake
  file(READ "${matio_SOURCE_DIR}/cmake/thirdParties.cmake" MATIO_THIRDPARTIES_CMAKE)
  
  # Find the MATIO_CREATE_ZLIB section
  string(FIND "${MATIO_THIRDPARTIES_CMAKE}" "macro(MATIO_CREATE_ZLIB" ZLIB_MACRO_POS)
  
  if(ZLIB_MACRO_POS EQUAL -1)
    message(FATAL_ERROR "Could not find MATIO_CREATE_ZLIB macro in thirdParties.cmake")
  endif()
  
  # Extract the header part before the MATIO_CREATE_ZLIB macro
  string(SUBSTRING "${MATIO_THIRDPARTIES_CMAKE}" 0 ${ZLIB_MACRO_POS} MATIO_THIRDPARTIES_HEADER)
  
  # Read our custom implementation
  configure_file("${CMAKE_CURRENT_SOURCE_DIR}/cmake/patch/matio_thirdParties.cmake.in" 
                "${CMAKE_CURRENT_BINARY_DIR}/matio_thirdParties.cmake.patched" @ONLY)
  file(READ "${CMAKE_CURRENT_BINARY_DIR}/matio_thirdParties.cmake.patched" MATIO_THIRDPARTIES_PATCHED)
  
  # Get everything after the if(MATIO_WITH_ZLIB) section
  string(FIND "${MATIO_THIRDPARTIES_CMAKE}" "endif()" MATIO_ZLIB_ENDIF_POS REVERSE)
  math(EXPR MATIO_FOOTER_POS "${MATIO_ZLIB_ENDIF_POS}+7") # 7 is the length of endif()
  string(SUBSTRING "${MATIO_THIRDPARTIES_CMAKE}" ${MATIO_FOOTER_POS} -1 MATIO_THIRDPARTIES_FOOTER)
  
  # Write the patched thirdParties.cmake
  file(WRITE "${matio_SOURCE_DIR}/cmake/thirdParties.cmake" 
       "${MATIO_THIRDPARTIES_HEADER}${MATIO_THIRDPARTIES_PATCHED}${MATIO_THIRDPARTIES_FOOTER}")
  
  message(STATUS "Successfully patched matio's thirdParties.cmake to use only static zlib")
endif()