#include <vector>
#include <cmath>
#include <Eigen/Dense>

// Retorna matriz S de rigidez local para elemento nodal linear (3x3)
std::vector<std::vector<double>> s_nodal(double x1, double y1,
                                         double x2, double y2,
                                         double x3, double y3) {
    std::vector<std::vector<double>> S(3, std::vector<double>(3, 0.0));

    double area = 0.5 * std::abs((x2 - x1)*(y3 - y1) - (x3 - x1)*(y2 - y1));

    // Monta matriz para obter os coeficientes da transformação baricêntrica
    Eigen::Matrix3d A;
    A << x1, x2, x3,
         y1, y2, y3,
         1.0, 1.0, 1.0;
    Eigen::Matrix3d Ainv = A.inverse();

    std::vector<double> b = { Ainv(1,0), Ainv(1,1), Ainv(1,2) };
    std::vector<double> c = { Ainv(2,0), Ainv(2,1), Ainv(2,2) };

    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            S[i][j] = area * (b[i]*b[j] + c[i]*c[j]);

    return S;
}
