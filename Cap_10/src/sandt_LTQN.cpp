#include <vector>
#include <cmath>
#include <utility>

// Dependências
std::vector<std::vector<double>> LTQN(const std::vector<double>& lambda,
                                      const std::vector<std::vector<double>>& nabla_lambda);

std::vector<double> curl_LTQN(const std::vector<double>& lambda,
                              const std::vector<std::vector<double>>& nabla_lambda);

// Quadratura de ordem 6: retorna pesos e pontos baricêntricos
std::pair<std::vector<double>, std::vector<std::vector<double>>> tri_quad(int order);

// Retorna S, T (8x8) para elemento LTQN com vértices (x1,y1), (x2,y2), (x3,y3)
std::pair<std::vector<std::vector<double>>, std::vector<std::vector<double>>>
sandt_LTQN(double x1, double y1, double x2, double y2, double x3, double y3) {
    std::vector<std::vector<double>> S(8, std::vector<double>(8, 0.0));
    std::vector<std::vector<double>> T(8, std::vector<double>(8, 0.0));

    // Área do triângulo
    double area = 0.5 * std::abs((x2 - x1)*(y3 - y1) - (x3 - x1)*(y2 - y1));

    // Cálculo de ∇λ_i (constantes no triângulo)
    Eigen::Matrix3d A;
    A << x1, x2, x3,
         y1, y2, y3,
         1.0, 1.0, 1.0;
    Eigen::Matrix3d Ainv = A.inverse();

    std::vector<std::vector<double>> nabla_lambda(3, std::vector<double>(2));
    for (int i = 0; i < 3; ++i) {
        nabla_lambda[i][0] = Ainv(1, i); // ∂λ/∂x
        nabla_lambda[i][1] = Ainv(2, i); // ∂λ/∂y
    }

    // Quadratura de ordem 6
    auto [weights, bary_coords] = tri_quad(6);
    int quad_pts = weights.size();

    for (int n = 0; n < quad_pts; ++n) {
        std::vector<double> lambda = bary_coords[n];
        std::vector<std::vector<double>> funcs = LTQN(lambda, nabla_lambda);
        std::vector<double> curls = curl_LTQN(lambda, nabla_lambda);

        for (int i = 0; i < 8; ++i) {
            for (int j = 0; j < 8; ++j) {
                double dot_f = funcs[i][0]*funcs[j][0] + funcs[i][1]*funcs[j][1];
                S[i][j] += weights[n] * curls[i] * curls[j];
                T[i][j] += weights[n] * dot_f;
            }
        }
    }

    // Multiplica pela área do triângulo
    for (int i = 0; i < 8; ++i)
        for (int j = 0; j < 8; ++j) {
            S[i][j] *= area;
            T[i][j] *= area;
        }

    return {S, T};
}
