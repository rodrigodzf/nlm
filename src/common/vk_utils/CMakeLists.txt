# Define common compile options
set(COMMON_RELEASE_OPTIONS -O3 -march=native -DNDEBUG -DEIGEN_USE_BLAS)

# Compiler-specific options
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    set(COMPILER_SPECIFIC_OPTIONS -ffast-math -finline-limit=1000000 -ffp-contract=fast)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    set(COMPILER_SPECIFIC_OPTIONS -ffast-math -ffp-contract=fast)
else()
    set(COMPILER_SPECIFIC_OPTIONS "")
endif()

# Combine options and wrap in generator expression for Release config
set(RELEASE_COMPILE_OPTIONS
    "$<$<CONFIG:Release>:${COMMON_RELEASE_OPTIONS};${COMPILER_SPECIFIC_OPTIONS}>"
)

set(LIB_NAME vk_utils)
# set(SRC_FILES 
    # vk_utils.cpp
# )

add_library(${LIB_NAME} INTERFACE)

# Apply the compile options to the interface library
target_compile_options(${LIB_NAME} INTERFACE ${RELEASE_COMPILE_OPTIONS})
target_include_directories(${LIB_NAME} INTERFACE include)
target_link_libraries(${LIB_NAME} INTERFACE Eigen3::Eigen BLAS::BLAS)