cmake_minimum_required(VERSION 3.16)

project(Cap_06 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(Eigen3 REQUIRED)

set(PROJECT_OUT_DIR "${CMAKE_SOURCE_DIR}/out")
file(MAKE_DIRECTORY "${PROJECT_OUT_DIR}")

add_library(cap06_common
    src/mom3d_globals.cpp
    src/trimesh3D.cpp
    src/edgemake_MoM.cpp
    src/outside_edge.cpp
    src/renumber_RWG.cpp
    src/edge_conx_elem.cpp
    src/find_local_dofs.cpp
    src/ComputeRho_c.cpp
    src/hat.cpp
    src/simplex_area.cpp
    src/tri_area3D.cpp
    src/GL_quad_rule.cpp
    src/tri_quad.cpp
    src/intg_sing_SGF.cpp
    src/Int_pq.cpp
    src/FillVVector.cpp
    src/FillZMatrixByEdge.cpp
    src/FillZMatrixByFace.cpp
    src/PostProcMoM.cpp
    src/sphereRCS.cpp
)

target_include_directories(cap06_common
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(cap06_common
    PUBLIC
        Eigen3::Eigen
)

target_compile_definitions(cap06_common
    PUBLIC
        PROJECT_OUT_DIR="${PROJECT_OUT_DIR}"
)

add_executable(MoM3D_demo MoM3D_demo.cpp)
target_link_libraries(MoM3D_demo PRIVATE cap06_common)

add_executable(converge_sphereRCS converge_sphereRCS.cpp)
target_link_libraries(converge_sphereRCS PRIVATE cap06_common)
