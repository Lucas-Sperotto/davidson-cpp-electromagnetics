#ifndef S_NODAL_H
#define S_NODAL_H

#include <vector>

// Retorna a matriz de rigidez local (3x3) para um elemento triangular de primeira ordem
std::vector<std::vector<double>> s_nodal(double x1, double y1,
                                         double x2, double y2,
                                         double x3, double y3);

#endif // S_NODAL_H
