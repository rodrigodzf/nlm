set(PACKAGE_NAME "modal_tilde")
set(PACKAGE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/${PACKAGE_NAME}")

# delete everything in the package directory
file(REMOVE_RECURSE ${PACKAGE_DIR})

# Copy the Max package
file(COPY
    ${CMAKE_CURRENT_SOURCE_DIR}/externals
    DESTINATION
    ${PACKAGE_DIR}
)

# Copy the extras directory
file(COPY
    ${CMAKE_CURRENT_SOURCE_DIR}/extras
    DESTINATION
    ${PACKAGE_DIR}
)

# Copy the package info, readme and license
file(COPY
    ${CMAKE_CURRENT_SOURCE_DIR}/package-info.json
    ${CMAKE_CURRENT_SOURCE_DIR}/README.md
    ${CMAKE_CURRENT_SOURCE_DIR}/LICENSE
    DESTINATION
    ${PACKAGE_DIR}
)

# Copy the help files
# file(COPY
#     "${CMAKE_CURRENT_SOURCE_DIR}/help/modal.plate~.maxhelp"
#     "${CMAKE_CURRENT_SOURCE_DIR}/help/modal.string~.maxhelp"
#     "${CMAKE_CURRENT_SOURCE_DIR}/help/mcs.modal.plate~.maxhelp"
#     "${CMAKE_CURRENT_SOURCE_DIR}/help/mcs.modal.string~.maxhelp"
#     DESTINATION
#     "${PACKAGE_DIR}/help"
# )

# Create symbolic links to help files
file(MAKE_DIRECTORY "${PACKAGE_DIR}/help")
execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink 
    "${CMAKE_CURRENT_SOURCE_DIR}/help/modal.plate~.maxhelp"
    "${PACKAGE_DIR}/help/modal.plate~.maxhelp")
execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink 
    "${CMAKE_CURRENT_SOURCE_DIR}/help/modal.string~.maxhelp"
    "${PACKAGE_DIR}/help/modal.string~.maxhelp")
execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink 
    "${CMAKE_CURRENT_SOURCE_DIR}/help/mcs.modal.plate~.maxhelp"
    "${PACKAGE_DIR}/help/mcs.modal.plate~.maxhelp")
execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink 
    "${CMAKE_CURRENT_SOURCE_DIR}/help/mcs.modal.string~.maxhelp"
    "${PACKAGE_DIR}/help/mcs.modal.string~.maxhelp")

# Copy the docs
file(COPY
    ${CMAKE_CURRENT_SOURCE_DIR}/docs
    DESTINATION
    ${PACKAGE_DIR}
)

# Copy the misc files
# file(COPY
#     ${CMAKE_CURRENT_SOURCE_DIR}/misc
#     DESTINATION
#     ${PACKAGE_DIR}
# )

# Create symbolic link to misc folder
execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink 
    ${CMAKE_CURRENT_SOURCE_DIR}/misc
    ${PACKAGE_DIR}/misc)

# execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink 
# ${CMAKE_CURRENT_SOURCE_DIR}/yaml
# ${PACKAGE_DIR}/yaml)