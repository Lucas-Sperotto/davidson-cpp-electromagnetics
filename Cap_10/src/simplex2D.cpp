#include <vector>
#include <cmath>
#include <cstdlib>

extern std::vector<std::vector<int>> ELEMENTS;            // ELEMENTS[e][0..2]
extern std::vector<std::vector<double>> NODE_COORD;       // NODE_COORD[n][0] = x, [1] = y

// Retorna (lambda, area) — coordenadas baricêntricas e área do triângulo
std::pair<std::vector<double>, double> simplex2D(int elem_num, double xp, double yp) {
    std::vector<int> trinodes = ELEMENTS[elem_num];

    double x1 = NODE_COORD[trinodes[0]][0], y1 = NODE_COORD[trinodes[0]][1];
    double x2 = NODE_COORD[trinodes[1]][0], y2 = NODE_COORD[trinodes[1]][1];
    double x3 = NODE_COORD[trinodes[2]][0], y3 = NODE_COORD[trinodes[2]][1];

    auto det = [](double a1, double a2, double a3,
                  double b1, double b2, double b3,
                  double c1, double c2, double c3) {
        return a1 * (b2 * c3 - b3 * c2)
             - a2 * (b1 * c3 - b3 * c1)
             + a3 * (b1 * c2 - b2 * c1);
    };

    double sigma  = det(1, x1, y1, 1, x2, y2, 1, x3, y3);
    double sigma1 = det(1, xp, yp, 1, x2, y2, 1, x3, y3);
    double sigma2 = det(1, x1, y1, 1, xp, yp, 1, x3, y3);
    double sigma3 = det(1, x1, y1, 1, x2, y2, 1, xp, yp);

    std::vector<double> lambda(3);
    lambda[0] = sigma1 / sigma;
    lambda[1] = sigma2 / sigma;
    lambda[2] = sigma3 / sigma;

    double tri_area = 0.5 * std::abs(sigma);

    return {lambda, tri_area};
}
