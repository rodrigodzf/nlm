# StaticLibraryConfig.cmake
# Centralized static library configuration for the entire project
#
# This module ensures all dependencies are built as static libraries
# to create self-contained Max externals.

# Force static libraries globally
set(BUILD_SHARED_LIBS OFF CACHE BOOL "Build shared libraries" FORCE)

# Additional static library preferences for common dependencies
set(MATIO_SHARED OFF CACHE BOOL "Build matio as shared library" FORCE)
set(ZLIB_ENABLE_SHARED OFF CACHE BOOL "Build zlib as shared library" FORCE)
set(HDF5_USE_STATIC_LIBRARIES ON CACHE BOOL "Use static HDF5 libraries" FORCE)
set(SKIP_INSTALL_SHARED_LIBRARIES ON CACHE BOOL "Skip installing shared libraries" FORCE)

# Ensure position independent code for static libraries
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

message(STATUS "Static library configuration: All libraries will be built statically")