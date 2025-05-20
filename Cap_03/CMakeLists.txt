cmake_minimum_required(VERSION 3.10)
project(Cap_03_FDTD_2D_3D)

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

include_directories(include)

# Executáveis
add_executable(fdtd_2d_demo src/fdtd_2d_demo.cpp
    src/gaussder.cpp
)

#target_compile_definitions(fdtd_2d_demo PRIVATE $<$<BOOL:${RUN_PML}>:RUN_PML>)

add_executable(fdtd3d
    main_fdtd3d.cpp
    src/fdtd3d.cpp
)

# Mensagem final
message(STATUS "")
message(STATUS "✅ Configuração concluída com sucesso!")
message(STATUS "")
message(STATUS "💡 Dica: compile com todos os núcleos disponíveis para maior desempenho:")
message(STATUS "   make -j$(nproc)")
message(STATUS "")