cmake_minimum_required(VERSION 3.16)

project(Cap_11 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(Eigen3 REQUIRED)

set(PROJECT_OUT_DIR "${CMAKE_SOURCE_DIR}/out")
file(MAKE_DIRECTORY "${PROJECT_OUT_DIR}")

add_library(cap11_common
    src/globals3d.cpp
    src/brick_mesh.cpp
    src/internal_tet_mesh.cpp
    src/read_gmsh2.cpp
    src/edgemake3D.cpp
    src/facemake3D.cpp
    src/free_dof3D.cpp
    src/free_dof3D_LTQN.cpp
    src/free_nodes3D.cpp
    src/renumber_dof3D.cpp
    src/renumber_dof3D_LTQN.cpp
    src/sandt3D.cpp
    src/sandt3D_LTQN.cpp
    src/tet_quad.cpp
    src/LTQN3D.cpp
    src/curl_LTQN3D.cpp
    src/TE_cavity_modes3D.cpp
    src/cavity_modes3D.cpp
    src/eig_err3D.cpp
    src/avg_mesh_length.cpp
)

target_include_directories(cap11_common
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(cap11_common
    PUBLIC
        Eigen3::Eigen
)

target_compile_definitions(cap11_common
    PUBLIC
        PROJECT_OUT_DIR="${PROJECT_OUT_DIR}"
)

add_executable(Eigen3D_v0 Eigen3D_v0.cpp)
target_link_libraries(Eigen3D_v0 PRIVATE cap11_common)

add_executable(Eigen3D_CTLN Eigen3D_CTLN.cpp)
target_link_libraries(Eigen3D_CTLN PRIVATE cap11_common)

add_executable(Eigen3D_LTQN Eigen3D_LTQN.cpp)
target_link_libraries(Eigen3D_LTQN PRIVATE cap11_common)

add_executable(FETD_FDTD FETD_FDTD.cpp)
target_link_libraries(FETD_FDTD PRIVATE cap11_common)

add_executable(test_tet_quad test_tet_quad.cpp)
target_link_libraries(test_tet_quad PRIVATE cap11_common)
