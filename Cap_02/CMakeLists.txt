cmake_minimum_required(VERSION 3.10)
project(Cap_02_FDTD_1D)

# Usa C++20
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Otimizações e flags específicas
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native -ffast-math")

# Diretórios
include_directories(${CMAKE_SOURCE_DIR}/src)
file(MAKE_DIRECTORY ${CMAKE_SOURCE_DIR}/out)

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

# Mensagem final
message(STATUS "")
message(STATUS "✅ Configuração concluída com sucesso!")
message(STATUS "")
message(STATUS "💡 Dica: compile com todos os núcleos disponíveis para maior desempenho:")
message(STATUS "   make -j$(nproc)")
message(STATUS "")