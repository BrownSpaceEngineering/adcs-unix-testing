cmake_minimum_required(VERSION 3.20)

project(ADCSTest)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)

set(CMAKE_C_FLAGS_DEBUG "-Wall -Wextra -Wpedantic -g" CACHE STRING "" FORCE)
set(CMAKE_C_FLAGS_RELEASE "-O3 -DNDEBUG" CACHE STRING "" FORCE)

# Sources shared by both executables
file(GLOB_RECURSE LINALG_SRC CONFIGURE_DEPENDS linalg/*.c)
set(COMMON_SRC
    ${LINALG_SRC}
    src/main.c
    src/test.c
    src/quat.c
    src/magnetosphere.c
    src/magnetosphere_update.c
)

set(COMMON_INCLUDES
    .
    linalg
    linalg/LinearAlgebra
    linalg/Lapack/Blas
    linalg/Lapack/F2c
    linalg/Lapack/Include
    linalg/Lapack/Install
    linalg/Lapack/Scr
)

# Original float implementation
add_executable(adcs-test ${COMMON_SRC})
target_include_directories(adcs-test PRIVATE ${COMMON_INCLUDES})
target_compile_definitions(adcs-test PRIVATE MAG_SPHERE_ORIGINAL)
target_link_libraries(adcs-test m)

# Precision-improved implementation (fixes 1-3)
add_executable(adcs-test-update ${COMMON_SRC})
target_include_directories(adcs-test-update PRIVATE ${COMMON_INCLUDES})
target_compile_definitions(adcs-test-update PRIVATE MAG_SPHERE_UPDATE)
target_link_libraries(adcs-test-update m)

# expose headers in an include directory if desired:
# target_include_directories(${EXECUTABLE} PUBLIC include)
