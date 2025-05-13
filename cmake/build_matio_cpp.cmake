# Set global shared library settings
set(BUILD_SHARED_LIBS OFF CACHE BOOL "Build shared libraries" FORCE)

# Define matio-cpp sources
set(MATIOCPP_SOURCE_DIR "${matio-cpp_SOURCE_DIR}")
set(MATIOCPP_SRC 
  "${MATIOCPP_SOURCE_DIR}/src/Variable.cpp"
  "${MATIOCPP_SOURCE_DIR}/src/ConversionUtilities.cpp"
  "${MATIOCPP_SOURCE_DIR}/src/MatvarHandler.cpp"
  "${MATIOCPP_SOURCE_DIR}/src/SharedMatvar.cpp"
  "${MATIOCPP_SOURCE_DIR}/src/WeakMatvar.cpp"
  "${MATIOCPP_SOURCE_DIR}/src/CellArray.cpp"
  "${MATIOCPP_SOURCE_DIR}/src/File.cpp"
  "${MATIOCPP_SOURCE_DIR}/src/Struct.cpp"
  "${MATIOCPP_SOURCE_DIR}/src/StructArray.cpp"
  "${MATIOCPP_SOURCE_DIR}/src/ExogenousConversions.cpp"
)

set(MATIOCPP_HDR 
  "${MATIOCPP_SOURCE_DIR}/include/matioCpp/Span.h"
  "${MATIOCPP_SOURCE_DIR}/include/matioCpp/VectorIterator.h"
  "${MATIOCPP_SOURCE_DIR}/include/matioCpp/ConversionUtilities.h"
  "${MATIOCPP_SOURCE_DIR}/include/matioCpp/ExogenousConversions.h"
  "${MATIOCPP_SOURCE_DIR}/include/matioCpp/EigenConversions.h"
  "${MATIOCPP_SOURCE_DIR}/include/matioCpp/Variable.h"
  "${MATIOCPP_SOURCE_DIR}/include/matioCpp/ForwardDeclarations.h"
  "${MATIOCPP_SOURCE_DIR}/include/matioCpp/Vector.h"
  "${MATIOCPP_SOURCE_DIR}/include/matioCpp/MultiDimensionalArray.h"
  "${MATIOCPP_SOURCE_DIR}/include/matioCpp/MatvarHandler.h"
  "${MATIOCPP_SOURCE_DIR}/include/matioCpp/SharedMatvar.h"
  "${MATIOCPP_SOURCE_DIR}/include/matioCpp/WeakMatvar.h"
  "${MATIOCPP_SOURCE_DIR}/include/matioCpp/Element.h"
  "${MATIOCPP_SOURCE_DIR}/include/matioCpp/CellArray.h"
  "${MATIOCPP_SOURCE_DIR}/include/matioCpp/File.h"
  "${MATIOCPP_SOURCE_DIR}/include/matioCpp/Struct.h"
  "${MATIOCPP_SOURCE_DIR}/include/matioCpp/StructArray.h"
  "${MATIOCPP_SOURCE_DIR}/include/matioCpp/StructArrayElement.h"
)

# Configure the header files
configure_file(
  "${MATIOCPP_SOURCE_DIR}/cmake/Config.h.in"
  "${CMAKE_BINARY_DIR}/matiocpp_build/include/matioCpp/Config.h"
  @ONLY
)

# Create the umbrella header
foreach(header ${MATIOCPP_HDR})
  get_filename_component(header_name ${header} NAME)
  set(include_command "#include <matioCpp/${header_name}>")
  set(umbrella_includes_list "${umbrella_includes_list}\n${include_command}")
endforeach()

configure_file("${MATIOCPP_SOURCE_DIR}/cmake/matioCpp.h.in"
               "${CMAKE_BINARY_DIR}/matiocpp_build/include/matioCpp/matioCpp.h" @ONLY)

list(APPEND MATIOCPP_HDR "${CMAKE_BINARY_DIR}/matiocpp_build/include/matioCpp/matioCpp.h")

# Add the matioCpp library manually
add_library(matioCpp STATIC ${MATIOCPP_SRC} ${MATIOCPP_HDR})

# Set include directories
target_include_directories(matioCpp
  PUBLIC
    "$<BUILD_INTERFACE:${MATIOCPP_SOURCE_DIR}/include>"
    "$<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/matiocpp_build/include>"
    "$<BUILD_INTERFACE:${matio_SOURCE_DIR}/src>"
    "$<BUILD_INTERFACE:${matio_configured_BINARY_DIR}/src>"
)

# Find the visit_struct dependency
include(FetchContent)
FetchContent_Declare(visit_struct
  GIT_REPOSITORY https://github.com/ami-iit/visit_struct
  GIT_TAG 47bc6a3aa7588a1f4db39579a0b6812569a76b56
)
FetchContent_MakeAvailable(visit_struct)

# Link dependencies
target_link_libraries(matioCpp 
  PUBLIC 
    MATIO::MATIO 
    visit_struct::visit_struct
)

if(Eigen3_FOUND)
  target_link_libraries(matioCpp PUBLIC Eigen3::Eigen)
  target_compile_definitions(matioCpp PUBLIC MATIOCPP_HAS_EIGEN)
endif()

# Set other properties
target_compile_features(matioCpp PUBLIC cxx_std_14)
set_target_properties(matioCpp PROPERTIES 
  VERSION 0.2.6
)

# Create an alias target
add_library(matioCpp::matioCpp ALIAS matioCpp)