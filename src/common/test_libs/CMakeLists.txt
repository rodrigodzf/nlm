set(EXE_NAME test_libs)
set(SRC_FILES main.cpp)

add_executable(${EXE_NAME} ${SRC_FILES})

# Set compiler flags
# target_compile_options(${EXE_NAME} PRIVATE
#     -Wall
#     -Wextra
#     -Werror
# )

# # Set include directories
# target_include_directories(${EXE_NAME} PRIVATE
#     ${CMAKE_CURRENT_SOURCE_DIR}/../include
# )

target_link_libraries(
    ${EXE_NAME}
    PUBLIC
    matioCpp::matioCpp
    cxxopts
    vk_utils
)
set_target_properties(${EXE_NAME} PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")
