#ifndef TRI_QUAD_H
#define TRI_QUAD_H

#include <vector>
#include <utility>

// Retorna pesos e coordenadas baricêntricas para quadratura de ordem 6 em triângulo
std::pair<std::vector<double>, std::vector<std::vector<double>>> tri_quad(int order);

#endif // TRI_QUAD_H
