#include <vector>
#include <utility>
#include <stdexcept>

// Retorna pesos e coordenadas baricêntricas para quadratura de ordem 6 em triângulo
std::pair<std::vector<double>, std::vector<std::vector<double>>> tri_quad(int order) {
    std::vector<double> w;
    std::vector<std::vector<double>> lambda;

    if (order == 6) {
        // Quadratura de Dunavant com 7 pontos (grau exato 5)
        // Fonte: Dunavant, D. A. (1985). High degree efficient symmetrical Gaussian quadrature rules
        const double a = 0.1012865073235;
        const double b = 0.4701420641051;
        const double w1 = 0.1125;
        const double w2 = 0.0661970763943;
        const double w3 = 0.0629695902724;

        // Ponto central
        w.push_back(w1);
        lambda.push_back({1.0/3, 1.0/3, 1.0/3});

        // 3 pontos simétricos com peso w2
        for (int i = 0; i < 3; ++i) {
            std::vector<double> lam(3, 0.0);
            lam[i] = b;
            lam[(i+1)%3] = 1.0 - 2.0*b;
            lam[(i+2)%3] = b;
            lambda.push_back(lam);
            w.push_back(w2);
        }

        // 3 pontos simétricos com peso w3
        for (int i = 0; i < 3; ++i) {
            std::vector<double> lam(3, 0.0);
            lam[i] = a;
            lam[(i+1)%3] = 1.0 - 2.0*a;
            lam[(i+2)%3] = a;
            lambda.push_back(lam);
            w.push_back(w3);
        }
    } else {
        throw std::runtime_error("Only order 6 quadrature implemented.");
    }

    return {w, lambda};
}
