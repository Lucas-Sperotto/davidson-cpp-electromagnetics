#include <vector>

// Tipo: vetor de 2 componentes (x, y)
using Vec2 = std::vector<double>;

// Retorna as 8 funções LTQN, cada uma sendo um vetor 2D
std::vector<Vec2> LTQN(const std::vector<double> &lambda, const std::vector<Vec2> &nabla_lambda)
{
    std::vector<Vec2> LTQN_functions(8, Vec2(2, 0.0));

    // Whitney 1-forms
    LTQN_functions[0][0] = lambda[0] * nabla_lambda[1][0] - lambda[1] * nabla_lambda[0][0];
    LTQN_functions[0][1] = lambda[0] * nabla_lambda[1][1] - lambda[1] * nabla_lambda[0][1];

    LTQN_functions[1][0] = lambda[0] * nabla_lambda[2][0] - lambda[2] * nabla_lambda[0][0];
    LTQN_functions[1][1] = lambda[0] * nabla_lambda[2][1] - lambda[2] * nabla_lambda[0][1];

    LTQN_functions[2][0] = lambda[1] * nabla_lambda[2][0] - lambda[2] * nabla_lambda[1][0];
    LTQN_functions[2][1] = lambda[1] * nabla_lambda[2][1] - lambda[2] * nabla_lambda[1][1];

    // LTLN functions (symmetric part)
    LTQN_functions[3][0] = lambda[0] * nabla_lambda[1][0] + lambda[1] * nabla_lambda[0][0];
    LTQN_functions[3][1] = lambda[0] * nabla_lambda[1][1] + lambda[1] * nabla_lambda[0][1];

    LTQN_functions[4][0] = lambda[0] * nabla_lambda[2][0] + lambda[2] * nabla_lambda[0][0];
    LTQN_functions[4][1] = lambda[0] * nabla_lambda[2][1] + lambda[2] * nabla_lambda[0][1];

    LTQN_functions[5][0] = lambda[1] * nabla_lambda[2][0] + lambda[2] * nabla_lambda[1][0];
    LTQN_functions[5][1] = lambda[1] * nabla_lambda[2][1] + lambda[2] * nabla_lambda[1][1];

    // Face functions
    LTQN_functions[6][0] = lambda[1] * lambda[2] * nabla_lambda[0][0] + lambda[0] * lambda[2] * nabla_lambda[1][0] - 2.0 * lambda[0] * lambda[1] * nabla_lambda[2][0];
    LTQN_functions[6][1] = lambda[1] * lambda[2] * nabla_lambda[0][1] + lambda[0] * lambda[2] * nabla_lambda[1][1] - 2.0 * lambda[0] * lambda[1] * nabla_lambda[2][1];

    LTQN_functions[7][0] = lambda[2] * lambda[0] * nabla_lambda[1][0] + lambda[1] * lambda[0] * nabla_lambda[2][0] - 2.0 * lambda[1] * lambda[2] * nabla_lambda[0][0];
    LTQN_functions[7][1] = lambda[2] * lambda[0] * nabla_lambda[1][1] + lambda[1] * lambda[0] * nabla_lambda[2][1] - 2.0 * lambda[1] * lambda[2] * nabla_lambda[0][1];

    return LTQN_functions;
}
