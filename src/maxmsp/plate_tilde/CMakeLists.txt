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
    matioCpp::matioCpp
    vk_utils
    Eigen3::Eigen
)

include(${C74_MIN_API_DIR}/script/min-posttarget.cmake)

#############################################################
# UNIT TEST
#############################################################

include(${C74_MIN_API_DIR}/test/min-object-unittest.cmake)

if(TARGET ${PROJECT_NAME}_test)

    target_include_directories(
        ${PROJECT_NAME}_test
        PUBLIC
        ${CMAKE_CURRENT_BINARY_DIR}
    )

    target_link_libraries(${PROJECT_NAME}_test
        PRIVATE
        Eigen3::Eigen
        matioCpp::matioCpp
        vk_utils
    )
endif()



# Copy libmatio from homebrew to the plate~ external's MacOS directory
if(APPLE)
    # Generate a script to properly set up the libraries
    set(INSTALL_NAME_SCRIPT "${CMAKE_BINARY_DIR}/fix_lib_paths.sh")
    file(WRITE ${INSTALL_NAME_SCRIPT} "#!/bin/bash
# Script to copy libmatio and fix install names

# The directory where the Max object resides
OBJECT_PATH=\"$1\"
MXOBJECT_PATH=\"$2\"
MATIO_DYLIB=\"${MATIO_DYLIB_PATH}\"
MATIO_FILENAME=$(basename \"$MATIO_DYLIB\")

# Create MacOS dir if needed
MACOS_DIR=\"${MXOBJECT_PATH}/Contents/MacOS\"
mkdir -p \"$MACOS_DIR\"

# Copy libmatio to the bundle
echo \"Copying $MATIO_DYLIB to $MACOS_DIR/\"
cp \"$MATIO_DYLIB\" \"$MACOS_DIR/\"

# Change libmatio's ID to use @loader_path
echo \"Setting libmatio ID\"
install_name_tool -id \"@loader_path/$MATIO_FILENAME\" \"$MACOS_DIR/$MATIO_FILENAME\"

# Find all dylibs in the Max object that reference libmatio
REFS=$(otool -L \"$OBJECT_PATH\" | grep -E \"libmatio.*dylib\" | awk '{print $1}')

# Change references in the Max object
for REF in $REFS; do
    echo \"Changing reference $REF to @loader_path/$MATIO_FILENAME in $OBJECT_PATH\"
    install_name_tool -change \"$REF\" \"@loader_path/$MATIO_FILENAME\" \"$OBJECT_PATH\"
done

# Check other dylibs in the bundle
for DYLIB in \"$MACOS_DIR\"/*.dylib; do
    if [ \"$DYLIB\" != \"$MACOS_DIR/$MATIO_FILENAME\" ]; then
        REFS=$(otool -L \"$DYLIB\" | grep -E \"libmatio.*dylib\" | awk '{print $1}')
        for REF in $REFS; do
            echo \"Changing reference $REF to @loader_path/$MATIO_FILENAME in $DYLIB\"
            install_name_tool -change \"$REF\" \"@loader_path/$MATIO_FILENAME\" \"$DYLIB\"
        done
    fi
done

echo \"Library paths fixed successfully.\"
")

    # Make the script executable
    file(CHMOD ${INSTALL_NAME_SCRIPT} PERMISSIONS 
         OWNER_READ OWNER_WRITE OWNER_EXECUTE 
         GROUP_READ GROUP_EXECUTE 
         WORLD_READ WORLD_EXECUTE)

    # Set the external bundle path
    set(PLATE_EXTERNAL_PATH "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${${PROJECT_NAME}_EXTERN_OUTPUT_NAME}.mxo")
    
    message(STATUS "PLATE_EXTERNAL_PATH: ${PLATE_EXTERNAL_PATH}")
    # Execute the fix script after build
    add_custom_command(
        TARGET ${PROJECT_NAME}
        POST_BUILD
        COMMAND "${INSTALL_NAME_SCRIPT}" "$<TARGET_FILE:${PROJECT_NAME}>" "${PLATE_EXTERNAL_PATH}"
        COMMENT "Setting up libmatio in plate~ external bundle"
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