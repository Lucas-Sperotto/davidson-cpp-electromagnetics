cmake_minimum_required(VERSION 3.10)
project(Cap_03_FDTD_2D_3D)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Diretórios
include_directories(${CMAKE_SOURCE_DIR}/src)
file(MAKE_DIRECTORY ${CMAKE_SOURCE_DIR}/Cap_03/out)

# Define macro para diretório de saída
add_compile_definitions(PROJECT_OUT_DIR="${CMAKE_SOURCE_DIR}/out")

include_directories(include)

# Executáveis
add_executable(fdtd_2d_demo src/fdtd_2d_demo.cpp
    src/gaussder.cpp
    src/fdtd2d_pml.cpp
)

target_compile_definitions(fdtd_2d_demo PRIVATE $<$<BOOL:${RUN_PML}>:RUN_PML>)

add_executable(fdtd3d
    main_fdtd3d.cpp
    src/fdtd3d.cpp
)
