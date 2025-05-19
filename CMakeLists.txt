cmake_minimum_required(VERSION 3.19)
project(modalsim VERSION 0.1.0)
message(STATUS "Version: ${PROJECT_VERSION}")

# We need to set a version for the package because the min-pretarget.cmake
# generates a new project per subdirectory
set(PACKAGE_VERSION ${PROJECT_VERSION})

# Extract package name from directory
string(REGEX REPLACE "(.*)/" "" THIS_PACKAGE_NAME "${CMAKE_CURRENT_SOURCE_DIR}")
message(STATUS "This package name is ${THIS_PACKAGE_NAME}")

# Enable testing support
enable_testing()

# C++17 standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_OSX_DEPLOYMENT_TARGET "10.15")
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Set static build preference globally
set(BUILD_SHARED_LIBS OFF CACHE BOOL "Build shared libraries" FORCE)

# macOS specific settings
if(APPLE)
    if(${CMAKE_GENERATOR} MATCHES "Xcode")
        if(${XCODE_VERSION} VERSION_LESS 10)
            message(STATUS "Xcode 10 or higher is required. Please install from the Mac Store.")
            return()
        elseif(${XCODE_VERSION} VERSION_GREATER_EQUAL 12)
            set(C74_BUILD_FAT YES)
        endif()
    endif()

    if(NOT CMAKE_OSX_ARCHITECTURES)
        if(C74_BUILD_FAT)
            set(CMAKE_OSX_ARCHITECTURES "x86_64;arm64" CACHE STRING "macOS architecture" FORCE)
        else()
            set(CMAKE_OSX_ARCHITECTURES ${CMAKE_SYSTEM_PROCESSOR} CACHE STRING "macOS architecture" FORCE)
        endif()

        message("CMAKE_OSX_ARCHITECTURES set to ${CMAKE_OSX_ARCHITECTURES}")
    endif()

    set(CMAKE_CXX_FLAGS "-faligned-allocation")

    if(CMAKE_OSX_ARCHITECTURES STREQUAL "")
        set(CMAKE_OSX_ARCHITECTURES ${CMAKE_HOST_SYSTEM_PROCESSOR})
    endif()

    message("CMAKE_OSX_ARCHITECTURES: ${CMAKE_OSX_ARCHITECTURES}")
endif()

#### Max SDK setup
set(C74_MIN_API_DIR ${CMAKE_CURRENT_SOURCE_DIR}/third_party/min-api)
include(${C74_MIN_API_DIR}/script/min-package.cmake)

# disable warnings for min-api
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    add_compile_options(
        -Wno-multichar 
        -Wno-non-template-friend
    )
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    add_compile_options(
        -Wno-deprecated-declarations
        -Wno-unused-private-field
    )
endif()

#### Max SDK setup

#### External dependencies

include(CPM.cmake)
set(LIB_DIR ${CMAKE_CURRENT_SOURCE_DIR}/third_party)

CPMAddPackage(
    NAME cxxopts
    GITHUB_REPOSITORY jarro2783/cxxopts
    VERSION 3.2.0
    OPTIONS
    "CXXOPTS_BUILD_EXAMPLES Off"
    "CXXOPTS_BUILD_TESTS Off"
)

# Find Eigen3 first, as it's needed by multiple components
find_package(Eigen3 REQUIRED)

# === zlib (static) ===
# First, just download zlib without configuring it
CPMAddPackage(
  NAME zlib
  GITHUB_REPOSITORY madler/zlib
  GIT_TAG v1.3.1
  DOWNLOAD_ONLY YES
)

# Create patched zlib CMakeLists.txt template
configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake/patch/zlib_CMakeLists.txt.in"
  "${CMAKE_CURRENT_BINARY_DIR}/zlib_CMakeLists.txt.in"
  @ONLY
)

# Patch the zlib CMakeLists.txt file before building
if(EXISTS "${zlib_SOURCE_DIR}/CMakeLists.txt")
  # Find the line where the library definitions start
  file(READ "${zlib_SOURCE_DIR}/CMakeLists.txt" ZLIB_ORIGINAL_CMAKE)
  string(FIND "${ZLIB_ORIGINAL_CMAKE}" "#============================================================================\n# zlib" ZLIB_START_POS)
  if(NOT ZLIB_START_POS EQUAL -1)
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
  endif()
endif()

# Now configure zlib with our patched CMakeLists.txt
CPMAddPackage(
  NAME zlib_configured
  SOURCE_DIR "${zlib_SOURCE_DIR}"
  OPTIONS
    "CMAKE_POSITION_INDEPENDENT_CODE ON"
    "ZLIB_BUILD_EXAMPLES OFF"
    "BUILD_SHARED_LIBS OFF"
    "SKIP_INSTALL_SHARED_LIBRARIES ON"
    "ZLIB_ENABLE_SHARED OFF"
)

# === matio (static) ===
set(MATIO_INSTALL_DIR ${CMAKE_BINARY_DIR}/_install/matio)
list(PREPEND CMAKE_PREFIX_PATH ${MATIO_INSTALL_DIR})

# Define policy to fix include directory issues
if(POLICY CMP0076)
  cmake_policy(SET CMP0076 NEW)
endif()

# Download matio first
CPMAddPackage(
  NAME matio
  GITHUB_REPOSITORY tbeu/matio
  GIT_TAG v1.5.28
  DOWNLOAD_ONLY YES
)

# Patch matio's thirdParties.cmake
if(EXISTS "${matio_SOURCE_DIR}/cmake/thirdParties.cmake")
  # Find the MATIO_CREATE_ZLIB macro in thirdParties.cmake
  file(READ "${matio_SOURCE_DIR}/cmake/thirdParties.cmake" MATIO_THIRDPARTIES_CMAKE)
  
  # Find the MATIO_CREATE_ZLIB section
  string(FIND "${MATIO_THIRDPARTIES_CMAKE}" "macro(MATIO_CREATE_ZLIB" ZLIB_MACRO_POS)
  
  if(NOT ZLIB_MACRO_POS EQUAL -1)
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
endif()

# Set static libraries in matio's CMakeLists.txt directly
if(EXISTS "${matio_SOURCE_DIR}/CMakeLists.txt")
  file(READ "${matio_SOURCE_DIR}/CMakeLists.txt" MATIO_CMAKE)
  # Add a line to force static libraries at the beginning of the file
  string(REPLACE "cmake_minimum_required" "set(BUILD_SHARED_LIBS OFF CACHE BOOL \"Build shared libraries\" FORCE)\ncmake_minimum_required" MATIO_CMAKE "${MATIO_CMAKE}")
  file(WRITE "${matio_SOURCE_DIR}/CMakeLists.txt" "${MATIO_CMAKE}")
  message(STATUS "Successfully patched matio's CMakeLists.txt to use only static libraries")
endif()

# Now configure matio with our patches
CPMAddPackage(
  NAME matio_configured
  SOURCE_DIR "${matio_SOURCE_DIR}"
  OPTIONS
    "CMAKE_POSITION_INDEPENDENT_CODE ON"
    "BUILD_SHARED_LIBS OFF"
    "MATIO_SHARED OFF"
    "MATIO_MAT73 OFF"
    "MATIO_WITH_HDF5 OFF"
    "MATIO_WITH_ZLIB ON"
    "MATIO_BUILD_TESTING OFF"
    "HDF5_USE_STATIC_LIBRARIES ON"
    "SKIP_INSTALL_SHARED_LIBRARIES ON"
    "ZLIB_INCLUDE_DIR=${zlib_SOURCE_DIR}"
    "ZLIB_LIBRARY=${CMAKE_BINARY_DIR}/_deps/zlib_configured-build/libz.a"
)

# Create the MATIO::MATIO alias target that matio-cpp expects
if (matio_configured_ADDED)
  if (NOT TARGET MATIO::MATIO)
    add_library(MATIO::MATIO ALIAS matio)
  endif()
  
  # Fix the include directories to use generator expressions 
  # This avoids the build directory prefixes that cause errors
  get_target_property(MATIO_INCLUDES matio INTERFACE_INCLUDE_DIRECTORIES)
  if(MATIO_INCLUDES)
    set_target_properties(matio PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "")
    target_include_directories(matio
      PUBLIC
        "$<BUILD_INTERFACE:${matio_SOURCE_DIR}/src>"
        "$<BUILD_INTERFACE:${matio_configured_BINARY_DIR}/src>"
    )
  endif()
endif()

# === matio-cpp (static) ===
# We need a fully custom build of matio-cpp to avoid installation issues
# First, download the source
CPMAddPackage(
  NAME matio-cpp
  GITHUB_REPOSITORY ami-iit/matio-cpp
  GIT_TAG v0.2.6
  VERSION 0.2.6
  OPTIONS
  DOWNLOAD_ONLY YES
)

include(cmake/build_matio_cpp.cmake)
# Add the cmake module path to find our custom modules
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules" ${CMAKE_MODULE_PATH})

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
find_package(ARPACK REQUIRED)

# Print status messages for dependencies
message(STATUS "Found EIGEN3: ${EIGEN3_INCLUDE_DIR}")
message(STATUS "Found BLAS: ${BLAS_LIBRARIES}")
message(STATUS "Found LAPACK: ${LAPACK_LIBRARIES}")
message(STATUS "Found ARPACK: ${ARPACK_LIBRARIES}")

# Set RPATH to look for libraries in the same directory as the executable
if(APPLE)
    set(CMAKE_MACOSX_RPATH TRUE)
    set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
    set(CMAKE_BUILD_RPATH_USE_ORIGIN TRUE)
    set(CMAKE_INSTALL_RPATH "@loader_path")
endif()

# Add subdirectories for your Max objects
add_subdirectory(src/maxmsp/modal.plate_tilde)
add_subdirectory(src/maxmsp/modal.string_tilde)
add_subdirectory(src/maxmsp/mcs.modal.plate_tilde)
add_subdirectory(src/maxmsp/mcs.modal.string_tilde)
add_subdirectory(src/common/vk_utils)
add_subdirectory(src/common/test_libs)
add_subdirectory(src/common/magpie)
add_subdirectory(src/common/coupling)

# Install the Max package
# include(cmake/make_package.cmake)