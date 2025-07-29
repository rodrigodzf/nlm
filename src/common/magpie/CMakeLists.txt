set(EXE_NAME test_magpie)
set(SRC_FILES main.cpp)

add_executable(${EXE_NAME} ${SRC_FILES})

target_include_directories(${EXE_NAME} PRIVATE include)

target_link_libraries(
    ${EXE_NAME}
    PUBLIC
    matioCpp::matioCpp
    cxxopts
    Eigen3::Eigen
    # LAPACK::LAPACK
    ${ARPACK_LIBRARIES}
)
set_target_properties(${EXE_NAME} PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")

# Define EIGEN_USE_ARPACK to enable ARPACK support
target_compile_definitions(${EXE_NAME} PRIVATE EIGEN_USE_ARPACK)

