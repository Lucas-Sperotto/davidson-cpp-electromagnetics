#ifndef CURL_LTQN_H
#define CURL_LTQN_H

#include <vector>

using Vec2 = std::vector<double>;

// Produto vetorial no plano (componente z)
double cross_z(const Vec2& a, const Vec2& b);

// Retorna os rotacionais das 8 funções LTQN (valores escalares em 2D)
std::vector<double> curl_LTQN(const std::vector<double>& lambda, const std::vector<Vec2>& nabla_lambda);

#endif // CURL_LTQN_H
