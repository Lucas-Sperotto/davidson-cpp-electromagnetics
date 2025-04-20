cmake_minimum_required(VERSION 3.10)
project(Cap_03_FDTD)

set(CMAKE_CXX_STANDARD 17)

include_directories(include)

add_executable(fdtd2d
    main.cpp
    src/fdtd2d.cpp
    src/gaussder.cpp
    src/fdtd2d_pml.cpp
)

target_compile_definitions(fdtd2d PRIVATE $<$<BOOL:${RUN_PML}>:RUN_PML>)

add_executable(fdtd3d
    src/fdtd3d.cpp
)
