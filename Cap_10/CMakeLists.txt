cmake_minimum_required(VERSION 3.10)
project(Capitulo10_CEM)

set(CMAKE_CXX_STANDARD 17)

# Caminhos
include_directories(include)
file(GLOB SOURCES "src/*.cpp")

# Executáveis principais
add_executable(Static2D Static2D.cpp ${SOURCES})
add_executable(Eigen2D Eigen2D.cpp ${SOURCES})
add_executable(Eigen2D_LTQN Eigen2D_LTQN.cpp ${SOURCES})

# Bibliotecas necessárias
find_package(Eigen3 REQUIRED)
target_link_libraries(Static2D Eigen3::Eigen)
target_link_libraries(Eigen2D Eigen3::Eigen)
target_link_libraries(Eigen2D_LTQN Eigen3::Eigen)
