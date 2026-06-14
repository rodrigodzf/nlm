set(PACKAGE_NAME "nlm")
set(PACKAGE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/${PACKAGE_NAME}")
set(PACKAGE_INFO_JSON "${CMAKE_CURRENT_BINARY_DIR}/package-info.json")

configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/package-info.json.in"
    "${PACKAGE_INFO_JSON}"
    @ONLY
)

# delete everything in the package directory
file(REMOVE_RECURSE ${PACKAGE_DIR})

# Option to choose between symbolic link and copy
option(USE_SYMLINKS "Use symbolic links instead of copying files" OFF)

# Function to handle file operations based on USE_SYMLINKS
function(handle_file_operation source dest)
    # Check if source exists
    if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${source}")
        if(USE_SYMLINKS)
            add_custom_target(link_${dest} ALL
                # First ensure parent directory exists
                COMMAND ${CMAKE_COMMAND} -E make_directory "${PACKAGE_DIR}"
                # Remove existing directory if it exists
                COMMAND ${CMAKE_COMMAND} -E remove_directory "${PACKAGE_DIR}/${dest}" || true
                # Create symbolic link
                COMMAND ${CMAKE_COMMAND} -E create_symlink
                    "${CMAKE_CURRENT_SOURCE_DIR}/${source}"
                    "${PACKAGE_DIR}/${dest}"
                COMMENT "Creating symbolic link for ${dest}"
            )
        else()
            add_custom_target(copy_${dest} ALL
                # First ensure parent directory exists
                COMMAND ${CMAKE_COMMAND} -E make_directory "${PACKAGE_DIR}"
                # Remove existing directory if it exists
                COMMAND ${CMAKE_COMMAND} -E remove_directory "${PACKAGE_DIR}/${dest}" || true
                # Copy directory
                COMMAND ${CMAKE_COMMAND} -E copy_directory
                    "${CMAKE_CURRENT_SOURCE_DIR}/${source}"
                    "${PACKAGE_DIR}/${dest}"
                COMMENT "Copying ${dest}"
            )
        endif()
    else()
        message(WARNING "Source directory ${source} does not exist, skipping...")
    endif()
endfunction()

# Handle all directories
handle_file_operation("externals" "externals")
handle_file_operation("extras" "extras")
handle_file_operation("docs" "docs")
handle_file_operation("misc" "misc")
handle_file_operation("help" "help")

if(TARGET copy_externals)
    foreach(EXTERNAL_TARGET
        nlm.plate_tilde
        nlm.string_tilde
        mcs.nlm.plate_tilde
        mcs.nlm.string_tilde
    )
        if(TARGET ${EXTERNAL_TARGET})
            add_dependencies(copy_externals ${EXTERNAL_TARGET})
        endif()
    endforeach()
endif()

# Handle individual files
if(USE_SYMLINKS)
    add_custom_target(link_individual_files ALL
        # First ensure parent directory exists
        COMMAND ${CMAKE_COMMAND} -E make_directory "${PACKAGE_DIR}"
        # Remove existing files if they exist
        COMMAND ${CMAKE_COMMAND} -E remove "${PACKAGE_DIR}/package-info.json" || true
        COMMAND ${CMAKE_COMMAND} -E remove "${PACKAGE_DIR}/README.md" || true
        COMMAND ${CMAKE_COMMAND} -E remove "${PACKAGE_DIR}/LICENSE" || true
        # Create symbolic links
        COMMAND ${CMAKE_COMMAND} -E create_symlink
            "${PACKAGE_INFO_JSON}"
            "${PACKAGE_DIR}/package-info.json"
        COMMAND ${CMAKE_COMMAND} -E create_symlink
            "${CMAKE_CURRENT_SOURCE_DIR}/README.md"
            "${PACKAGE_DIR}/README.md"
        COMMAND ${CMAKE_COMMAND} -E create_symlink
            "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE"
            "${PACKAGE_DIR}/LICENSE"
        COMMENT "Creating symbolic links for individual files"
    )
else()
    add_custom_target(copy_individual_files ALL
        # First ensure parent directory exists
        COMMAND ${CMAKE_COMMAND} -E make_directory "${PACKAGE_DIR}"
        # Copy files
        COMMAND ${CMAKE_COMMAND} -E copy
            "${PACKAGE_INFO_JSON}"
            "${CMAKE_CURRENT_SOURCE_DIR}/README.md"
            "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE"
            "${PACKAGE_DIR}"
        COMMENT "Copying individual files"
    )
endif()
