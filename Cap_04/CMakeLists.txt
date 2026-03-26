cmake_minimum_required(VERSION 3.12)
project(Cap_04_MoM)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Eigen3 3.3 REQUIRED NO_MODULE)

include_directories(${CMAKE_SOURCE_DIR}/include)
file(MAKE_DIRECTORY ${CMAKE_SOURCE_DIR}/out)
add_compile_definitions(PROJECT_OUT_DIR="${CMAKE_SOURCE_DIR}/out")

add_library(cap04_common STATIC
    src/mom_tm_solver.cpp
    src/cyl_tm_echo_width.cpp
)
target_link_libraries(cap04_common PRIVATE Eigen3::Eigen)

add_executable(MoM_2D_TM MoM_2D_TM.cpp)
target_link_libraries(MoM_2D_TM PRIVATE cap04_common Eigen3::Eigen)

add_executable(static_mom static_mom.cpp)
target_link_libraries(static_mom PRIVATE Eigen3::Eigen)

add_executable(thin_dipole thin_dipole.cpp)
target_link_libraries(thin_dipole PRIVATE Eigen3::Eigen)
