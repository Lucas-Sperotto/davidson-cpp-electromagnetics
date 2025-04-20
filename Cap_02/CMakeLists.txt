cmake_minimum_required(VERSION 3.10)
project(Cap_02_FDTD_1D)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Diretórios
include_directories(${CMAKE_SOURCE_DIR}/src)
file(MAKE_DIRECTORY ${CMAKE_SOURCE_DIR}/Cap_02/out)

# Define macro para diretório de saída
add_compile_definitions(PROJECT_OUT_DIR="${CMAKE_SOURCE_DIR}/out")

# FFTW
find_package(PkgConfig REQUIRED)
pkg_check_modules(FFTW REQUIRED fftw3)
include_directories(${FFTW_INCLUDE_DIRS})
link_directories(${FFTW_LIBRARY_DIRS})

# Executáveis
add_executable(fdtd_1D_demo src/fdtd_1D_demo.cpp)
target_link_libraries(fdtd_1D_demo ${FFTW_LIBRARIES} m)

add_executable(fdtd_1D_WB_demo src/fdtd_1D_WB_demo.cpp)
target_link_libraries(fdtd_1D_WB_demo ${FFTW_LIBRARIES} m)