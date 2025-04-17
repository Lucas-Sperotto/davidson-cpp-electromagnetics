#include <vector>

// Tipo: vetor 2D (x, y)
using Vec2 = std::vector<double>;

// Cross product (em 2D, com vetor z implícito)
double cross_z(const Vec2& a, const Vec2& b) {
    return a[0]*b[1] - a[1]*b[0];
}

// Retorna os rotacionais (escalares no plano 2D) das 8 funções LTQN
std::vector<double> curl_LTQN(const std::vector<double>& lambda, const std::vector<Vec2>& nabla_lambda) {
    std::vector<double> curl(8, 0.0);

    // Whitney 1-forms
    curl[0] = 2.0 * cross_z(nabla_lambda[0], nabla_lambda[1]);
    curl[1] = 2.0 * cross_z(nabla_lambda[0], nabla_lambda[2]);
    curl[2] = 2.0 * cross_z(nabla_lambda[1], nabla_lambda[2]);

    // LTLN (symmetric) functions: curl = 0
    curl[3] = 0.0;
    curl[4] = 0.0;
    curl[5] = 0.0;

    // Face functions
    curl[6] = -3.0 * lambda[1] * cross_z(nabla_lambda[0], nabla_lambda[2])
              -3.0 * lambda[0] * cross_z(nabla_lambda[1], nabla_lambda[2]);

    curl[7] = -3.0 * lambda[2] * cross_z(nabla_lambda[1], nabla_lambda[0])
              -3.0 * lambda[1] * cross_z(nabla_lambda[2], nabla_lambda[0]);

    return curl;
}
