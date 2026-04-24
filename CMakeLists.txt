cmake_minimum_required(VERSION 3.20)

project(ADCSTest)
set(EXECUTABLE adcs-test)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)

# glob all of linalg
file(GLOB_RECURSE SRC CONFIGURE_DEPENDS linalg/*.c)
# manual specified sources:
list(APPEND SRC src/main.c src/test.c src/quat.c src/laextension.c src/bdot.c src/body.c src/down_quat.c src/ecef2eci.c src/iterate.c src/kepler.c src/laextension.c src/magnetosphere.c src/moments2currents.c src/orbital2eci.c src/pd.c src/photodiode_determination.c src/quest.c src/sunvec.c src/torque2moments.c)
set(CMAKE_C_FLAGS_DEBUG "-Wall -Wextra -Wpedantic -g" CACHE STRING "" FORCE)
set(CMAKE_C_FLAGS_RELEASE "-O3 -DNDEBUG" CACHE STRING "" FORCE)

add_executable(${EXECUTABLE} ${SRC})
target_include_directories(${EXECUTABLE} PRIVATE .)
target_include_directories(${EXECUTABLE} PRIVATE linalg
    linalg/LinearAlgebra
    linalg/Lapack/Blas
    linalg/Lapack/F2c
    linalg/Lapack/Include
    linalg/Lapack/Install
    linalg/Lapack/Scr
)

target_link_libraries(adcs-test m)

# expose headers in an include directory if desired:
# target_include_directories(${EXECUTABLE} PUBLIC include)
