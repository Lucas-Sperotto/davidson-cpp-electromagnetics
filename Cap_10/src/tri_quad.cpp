#include <stdexcept>
#include <utility>
#include <vector>

// Retorna pesos e coordenadas baricentricas seguindo exatamente tri_quad.m.
std::pair<std::vector<double>, std::vector<std::vector<double>>> tri_quad(int order)
{
    if (order != 6)
        throw std::runtime_error("Only order 6 quadrature implemented.");

    std::vector<double> w(6, 0.0);
    std::vector<double> alpha(6, 0.0);
    std::vector<double> beta(6, 0.0);
    std::vector<double> gamma(6, 0.0);

    w[0] = 0.223381589678011;
    w[1] = w[0];
    w[2] = w[0];
    w[3] = 0.109951743655322;
    w[4] = w[3];
    w[5] = w[3];

    alpha[0] = 0.108103018168070;
    beta[0] = 0.445948490915965;
    gamma[0] = beta[0];

    alpha[1] = beta[0];
    alpha[2] = beta[0];
    beta[1] = alpha[0];
    beta[2] = beta[0];
    gamma[1] = beta[0];
    gamma[2] = alpha[0];

    alpha[3] = 0.816847572980459;
    beta[3] = 0.091576213509771;
    gamma[3] = beta[3];

    alpha[4] = beta[3];
    alpha[5] = beta[3];
    beta[4] = alpha[3];
    beta[5] = beta[3];
    gamma[4] = beta[3];
    gamma[5] = alpha[3];

    std::vector<std::vector<double>> lambda(6, std::vector<double>(3, 0.0));
    for (int i = 0; i < 6; ++i)
    {
        lambda[i][0] = alpha[i];
        lambda[i][1] = beta[i];
        lambda[i][2] = gamma[i];
    }

    return {w, lambda};
}
