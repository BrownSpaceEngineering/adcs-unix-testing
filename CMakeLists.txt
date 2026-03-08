cmake_minimum_required(VERSION 3.20)

project(ADCSTest)
set(EXECUTABLE adcs-test)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)

# to glob all src files in src/
file(GLOB_RECURSE SRC CONFIGURE_DEPENDS linalg/*.c)
# or manually specify:
set(SRC SRC test.c)
set(CMAKE_C_FLAGS_DEBUG "-Wall -Wextra -Wpedantic -g" CACHE STRING "" FORCE)
set(CMAKE_C_FLAGS_RELEASE "-O3 -DNDEBUG" CACHE STRING "" FORCE)

add_executable(${EXECUTABLE} ${SRC})
target_include_directories(${EXECUTABLE} PRIVATE .)
target_include_directories(${EXECUTABLE} PRIVATE linalg)
target_include_directories(${EXECUTABLE} PRIVATE linalg/LinearAlgebra)
target_include_directories(${EXECUTABLE} PRIVATE linalg/Lapack/Blas)
target_include_directories(${EXECUTABLE} PRIVATE linalg/Lapack/F2c)
target_include_directories(${EXECUTABLE} PRIVATE linalg/Lapack/Include)
target_include_directories(${EXECUTABLE} PRIVATE linalg/Lapack/Install)
target_include_directories(${EXECUTABLE} PRIVATE linalg/Lapack/Scr)
# TODO: expose headers in an include directory if desired:
# target_include_directories(${EXECUTABLE} PUBLIC include)
