#include <vector>

#include "simplex2D.h"

// Retorna as funções de base nodais (alpha_i) no ponto (x, y) dentro do elemento
std::vector<double> nodal_1st(int elem_num, double x, double y) {
    std::vector<double> lambda = simplex2D(elem_num, x, y).first; // coordenadas baricêntricas
    return lambda; // as funções nodais são diretamente as coordenadas de baricentro
}
