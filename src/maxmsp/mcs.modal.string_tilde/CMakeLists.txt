include(${C74_MIN_API_DIR}/script/min-pretarget.cmake)

include_directories( 
	"${C74_INCLUDES}"
)

set( SOURCE_FILES
	${PROJECT_NAME}.cpp
)

add_library( 
	${PROJECT_NAME} 
	MODULE
	${SOURCE_FILES}
)

# Comment this out as it will be automatically set by the directory name
# set(${PROJECT_NAME}_EXTERN_OUTPUT_NAME "mcs.modal.string~")

message(STATUS "Version in CMakeLists.txt: ${PACKAGE_VERSION}")

# Generate version header file
configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/version.h.in
    ${CMAKE_CURRENT_BINARY_DIR}/version.h
    @ONLY
)

# Add binary dir to include path to find the generated header
target_include_directories(
    ${PROJECT_NAME}
    PUBLIC
    ${CMAKE_CURRENT_BINARY_DIR}
)

target_compile_definitions(
    ${PROJECT_NAME}
    PUBLIC
    HAS_VERSION_H
)

# Set up RPATH to find dylibs in the same directory
# if(APPLE)
#     set_target_properties(${PROJECT_NAME} PROPERTIES
#         INSTALL_RPATH "@loader_path"
#         BUILD_WITH_INSTALL_RPATH TRUE
#     )
# endif()

# Forcibly set RPATH on the target
if(APPLE)
    set_target_properties(${PROJECT_NAME} PROPERTIES
        BUILD_WITH_INSTALL_RPATH TRUE
        INSTALL_RPATH "@loader_path"
        INSTALL_NAME_DIR "@loader_path"
    )
endif()

# Link libraries and their include directories will be propagated automatically
target_link_libraries(
    ${PROJECT_NAME}
    PUBLIC
    vk_utils
    Eigen3::Eigen
)

include(${C74_MIN_API_DIR}/script/min-posttarget.cmake)

#############################################################
# UNIT TEST
#############################################################

# include(${C74_MIN_API_DIR}/test/min-object-unittest.cmake)

# if(TARGET ${PROJECT_NAME}_test)

#     target_include_directories(
#         ${PROJECT_NAME}_test
#         PUBLIC
#         ${CMAKE_CURRENT_BINARY_DIR}
#     )

#     target_link_libraries(${PROJECT_NAME}_test
#         PRIVATE
#         Eigen3::Eigen
#         vk_utils
#     )
# endif()

# Copy libmatio from homebrew to the plate~ external's MacOS directory
if(APPLE)
    # Set the external bundle path
    set(PLATE_EXTERNAL_PATH "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${${PROJECT_NAME}_EXTERN_OUTPUT_NAME}.mxo")
    
    message(STATUS "PLATE_EXTERNAL_PATH: ${PLATE_EXTERNAL_PATH}")
    # Execute the fix script after build
    add_custom_command(
        TARGET ${PROJECT_NAME}
        POST_BUILD
        COMMAND "${INSTALL_NAME_SCRIPT}" "$<TARGET_FILE:${PROJECT_NAME}>" "${PLATE_EXTERNAL_PATH}"
        COMMENT "Setting up libmatio in mcs.modal.string~ external bundle"
    )

    if (CMAKE_OSX_ARCHITECTURES STREQUAL "arm64")
        add_custom_command( 
            TARGET ${PROJECT_NAME} 
            POST_BUILD 
            COMMAND "codesign" "--force" "--deep" "-s" "-" "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${${PROJECT_NAME}_EXTERN_OUTPUT_NAME}.mxo"
            COMMENT "Codesign external" 
         )
    endif()
endif()