cmake_minimum_required(VERSION 3.12)
project(Cap_09_FEM_1D)

# Usa C++20
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Otimizações e flags específicas
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native")

find_package(Eigen3 3.3 REQUIRED NO_MODULE)

add_executable(fem_tl
    FEM_1D_1st.cpp
    FEM_1D_solver.cpp
    FEM_pp.cpp
)

target_link_libraries(fem_tl PRIVATE Eigen3::Eigen)

# Mensagem final
message(STATUS "")
message(STATUS "✅ Configuração concluída com sucesso!")
message(STATUS "")
message(STATUS "💡 Dica: compile com todos os núcleos disponíveis para maior desempenho:")
message(STATUS "   make -j$(nproc)")
message(STATUS "")