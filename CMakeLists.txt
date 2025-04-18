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

# macOS specific settings
if(APPLE)
    if(${CMAKE_GENERATOR} MATCHES "Xcode")
        if(${XCODE_VERSION} VERSION_LESS 10)
            message(STATUS "Xcode 10 or higher is required. Please install from the Mac App Store.")
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

CPMAddPackage(
    NAME matio-cpp
    GIT_TAG v0.2.6
    VERSION 0.2.6
    GITHUB_REPOSITORY ami-iit/matio-cpp
    SOURCE_DIR ${LIB_DIR}/matio-cpp
)

find_package(Eigen3 REQUIRED)
find_package(BLAS REQUIRED)
# find_package(LAPACK REQUIRED)
# find_package(ARPACK REQUIRED)
#### External dependencies

# Add the cmake module path to find our custom modules
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules" ${CMAKE_MODULE_PATH})

# Find ARPACK using our custom module
find_package(LAPACK REQUIRED)
find_package(ARPACK REQUIRED)

# Print status messages for dependencies
message(STATUS "Found EIGEN3: ${EIGEN3_INCLUDE_DIR}")
message(STATUS "Found BLAS: ${BLAS_LIBRARIES}")
message(STATUS "Found LAPACK: ${LAPACK_LIBRARIES}")
message(STATUS "Found ARPACK: ${ARPACK_LIBRARIES}")

# Set RPATH
# set(CMAKE_MACOSX_RPATH TRUE)
# set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
# set(CMAKE_BUILD_RPATH_USE_ORIGIN TRUE)
# set(CMAKE_INSTALL_RPATH "@loader_path")

# Add subdirectories for your Max objects
add_subdirectory(src/maxmsp/plate_tilde)
add_subdirectory(src/common/vk_utils)
add_subdirectory(src/common/test_libs)
# add_subdirectory(src/common/magpie)
# Install configuration
# Install the Max package
install(
    DIRECTORY
    ${CMAKE_CURRENT_SOURCE_DIR}/externals
    DESTINATION
    ${CMAKE_CURRENT_SOURCE_DIR}/${THIS_PACKAGE_NAME}
)

# install the extras directory
install(
    DIRECTORY
    ${CMAKE_CURRENT_SOURCE_DIR}/extras
    DESTINATION
    ${CMAKE_CURRENT_SOURCE_DIR}/${THIS_PACKAGE_NAME}
)

# install the package info, readme and license
install(
    FILES
    ${CMAKE_CURRENT_SOURCE_DIR}/package-info.json
    ${CMAKE_CURRENT_SOURCE_DIR}/README.md
    ${CMAKE_CURRENT_SOURCE_DIR}/LICENSE
    DESTINATION
    ${CMAKE_CURRENT_SOURCE_DIR}/${THIS_PACKAGE_NAME}
)

# Install the help files (uncomment when you have help files)
# install(
# FILES
# ${CMAKE_CURRENT_SOURCE_DIR}/src/maxmsp/modal_resonator/modal_resonator.maxhelp
# DESTINATION
# ${CMAKE_CURRENT_SOURCE_DIR}/${THIS_PACKAGE_NAME}/help
# )