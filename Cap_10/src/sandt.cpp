#include <vector>
#include <cmath>
#include <array>
#include <cassert>
#include <Eigen/Dense>

// Estrutura auxiliar: tripla de índices das arestas locais
extern std::vector<std::array<int, 2>> LOCALEDGENODES; // geralmente = {{0,1},{0,2},{1,2}}

// Retorna S e T: matrizes locais 3×3 do elemento de Whitney
std::pair<std::vector<std::vector<double>>, std::vector<std::vector<double>>>
sandt(double x1, double y1, double x2, double y2, double x3, double y3) {
    std::vector<std::vector<double>> S(3, std::vector<double>(3, 0.0));
    std::vector<std::vector<double>> T(3, std::vector<double>(3, 0.0));

    double area = 0.5 * std::abs((x2 - x1)*(y3 - y1) - (x3 - x1)*(y2 - y1));

    // Matriz de transformação para coordenadas baricêntricas
    Eigen::Matrix3d A;
    A << x1, x2, x3,
         y1, y2, y3,
         1.0, 1.0, 1.0;
    Eigen::Matrix3d Ainv = A.inverse();

    std::vector<double> a = { Ainv(0,0), Ainv(0,1), Ainv(0,2) };
    std::vector<double> b = { Ainv(1,0), Ainv(1,1), Ainv(1,2) };
    std::vector<double> c = { Ainv(2,0), Ainv(2,1), Ainv(2,2) };

    std::vector<std::vector<double>> phi(3, std::vector<double>(3));
    std::vector<std::vector<double>> v_z(3, std::vector<double>(3));

    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) {
            phi[i][j] = b[i]*b[j] + c[i]*c[j];
            v_z[i][j] = b[i]*c[j] - b[j]*c[i];
        }

    const double M[3][3] = {
        {2.0/12, 1.0/12, 1.0/12},
        {1.0/12, 2.0/12, 1.0/12},
        {1.0/12, 1.0/12, 2.0/12}
    };

    for (int ii = 0; ii < 3; ++ii) {
        for (int jj = 0; jj < 3; ++jj) {
            int i1 = LOCALEDGENODES[ii][0];
            int i2 = LOCALEDGENODES[ii][1];
            int j1 = LOCALEDGENODES[jj][0];
            int j2 = LOCALEDGENODES[jj][1];

            S[ii][jj] = 4 * area * v_z[i1][i2] * v_z[j1][j2];
            T[ii][jj] = area * ( phi[i2][j2]*M[i1][j1] - phi[i2][j1]*M[i1][j2]
                               - phi[i1][j2]*M[i2][j1] + phi[i1][j1]*M[i2][j2] );
        }
    }

    return {S, T};
}
