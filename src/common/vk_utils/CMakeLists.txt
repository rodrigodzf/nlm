# Define the compile options
# set(
    # RELEASE_COMPILE_OPTIONS
    # "$<$<CONFIG:Release>:-O3;-mfma;-march=native;-DNDEBUG;-DEIGEN_USE_BLAS;-ffast-math;-finline-limit=1000000;-ffp-contract=fast>"
# )

set(LIB_NAME vk_utils)
# set(SRC_FILES 
    # vk_utils.cpp
# )

add_library(${LIB_NAME} INTERFACE)

# Apply the compile options to the interface library
target_compile_options(${LIB_NAME} INTERFACE ${RELEASE_COMPILE_OPTIONS})
target_include_directories(${LIB_NAME} INTERFACE include)
target_link_libraries(${LIB_NAME} INTERFACE Eigen3::Eigen)