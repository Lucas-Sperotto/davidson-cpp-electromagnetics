#include <vector>
#include <Eigen/Dense>

// Externas
extern std::vector<std::vector<int>> ELEMENTS;
extern std::vector<std::vector<double>> NODE_COORD;

// Retorna as 3 funções de Whitney no triângulo elem_num, avaliadas no ponto (x,y)
std::vector<std::vector<double>> whitney(int elem_num, double x, double y) {
    const auto& trinodes = ELEMENTS[elem_num];

    double x1 = NODE_COORD[trinodes[0]][0], y1 = NODE_COORD[trinodes[0]][1];
    double x2 = NODE_COORD[trinodes[1]][0], y2 = NODE_COORD[trinodes[1]][1];
    double x3 = NODE_COORD[trinodes[2]][0], y3 = NODE_COORD[trinodes[2]][1];

    // Calcula coordenadas baricêntricas
    auto [lambda, _] = simplex2D(elem_num, x, y);

    // Matriz para ∇λ_i
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

    std::vector<std::vector<double>> w(3, std::vector<double>(2, 0.0));

    // Whitney 1-forms (associação hardcoded: edge 12, 13, 23)
    for (int i = 0; i < 2; ++i) {
        w[0][i] = lambda[0]*nabla_lambda[1][i] - lambda[1]*nabla_lambda[0][i]; // edge 1–2
        w[1][i] = lambda[0]*nabla_lambda[2][i] - lambda[2]*nabla_lambda[0][i]; // edge 1–3
        w[2][i] = lambda[1]*nabla_lambda[2][i] - lambda[2]*nabla_lambda[1][i]; // edge 2–3
    }

    return w; // 3 funções de base vetoriais (2 componentes cada)
}
