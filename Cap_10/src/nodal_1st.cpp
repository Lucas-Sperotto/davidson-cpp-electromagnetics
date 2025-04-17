#include <vector>

// Assinatura da função simplex2D (deve já estar implementada)
std::vector<double> simplex2D(int elem_num, double x, double y);

// Retorna as funções de base nodais (alpha_i) no ponto (x, y) dentro do elemento
std::vector<double> nodal_1st(int elem_num, double x, double y) {
    std::vector<double> lambda = simplex2D(elem_num, x, y);
    return lambda; // as funções nodais são diretamente as coordenadas de baricentro
}
