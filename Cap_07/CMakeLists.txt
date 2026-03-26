cmake_minimum_required(VERSION 3.16)

project(Cap_07 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(Eigen3 REQUIRED)

set(PROJECT_OUT_DIR "${CMAKE_SOURCE_DIR}/out")
file(MAKE_DIRECTORY "${PROJECT_OUT_DIR}")

add_library(cap07_common
    src/chapter7_common.cpp
    src/D_TE.cpp
    src/D_TM.cpp
    src/F.cpp
    src/F_reg1.cpp
    src/F_reg2.cpp
    src/F_reg3.cpp
    src/F_static.cpp
    src/root_D_TM.cpp
    src/V_int.cpp
)

target_include_directories(cap07_common
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(cap07_common
    PUBLIC
        Eigen3::Eigen
)

target_compile_definitions(cap07_common
    PUBLIC
        PROJECT_OUT_DIR="${PROJECT_OUT_DIR}"
)

add_executable(scalar_pot scalar_pot.cpp)
target_link_libraries(scalar_pot PRIVATE cap07_common)

add_executable(V_pot_eps V_pot_eps.cpp)
target_link_libraries(V_pot_eps PRIVATE cap07_common)

add_executable(V_pot_height V_pot_height.cpp)
target_link_libraries(V_pot_height PRIVATE cap07_common)

add_executable(MoM_Som MoM_Som.cpp)
target_link_libraries(MoM_Som PRIVATE cap07_common)
