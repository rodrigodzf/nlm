include(${C74_MIN_API_DIR}/script/min-pretarget.cmake)

# Add the visit_struct paths - this is found in the build directory
if(NOT TARGET visit_struct::visit_struct)
    find_package(visit_struct QUIET)
    if(NOT visit_struct_FOUND)
        message(STATUS "Using visit_struct from build directory")
        include_directories(${CMAKE_BINARY_DIR}/_deps/visit_struct-src/include)
    endif()
endif()

# Find Eigen
# find_package(Eigen3 REQUIRED)

# # Name of your external
# set(C74_LIBRARY_NAME modal_resonator)

# # Source files for your external
# set(PROJECT_SRC
#     "modal_resonator.cpp"
# )

# # Header files for your external
# set(PROJECT_HEADERS
#     "modal_resonator.h"
# )

# Only keep the C74_INCLUDES here which is specific to Max/MSP
include_directories( 
	"${C74_INCLUDES}"
)

set( SOURCE_FILES
	${PROJECT_NAME}.cpp
)

# Define the main Max/MSP external
add_library( 
	${PROJECT_NAME} 
	MODULE
	${SOURCE_FILES}
)

target_compile_definitions(
    ${PROJECT_NAME}
    PUBLIC
    VERSION="${CMAKE_PROJECT_VERSION}"
)

# Link libraries and their include directories will be propagated automatically
target_link_libraries(
    ${PROJECT_NAME}
    PUBLIC
    matioCpp::matioCpp
    vk_utils
    Eigen3::Eigen
)

# Only include the vk_utils/include directory explicitly if needed
# target_include_directories(${PROJECT_NAME} 
#     PUBLIC
#     ${CMAKE_SOURCE_DIR}/src/common/vk_utils/include
# )

include(${C74_MIN_API_DIR}/script/min-posttarget.cmake)

#############################################################
# UNIT TEST
#############################################################

include(${C74_MIN_API_DIR}/test/min-object-unittest.cmake)

# The test might not inherit the include directories properly, so set them explicitly
if(TARGET ${PROJECT_NAME}_test)
    target_include_directories(${PROJECT_NAME}_test 
        PRIVATE
        ${EIGEN3_INCLUDE_DIR}
        ${CMAKE_SOURCE_DIR}/src/common/vk_utils/include
        ${CMAKE_BINARY_DIR}/_deps/visit_struct-src/include
    )
    
    target_link_libraries(${PROJECT_NAME}_test
        PRIVATE
        Eigen3::Eigen
        matioCpp::matioCpp
        vk_utils
    )
endif()

# # Set compiler options for the test
# target_compile_options(
#     "${PROJECT_NAME}_test"
#     PUBLIC
#     -Wno-deprecated-declarations
# )
