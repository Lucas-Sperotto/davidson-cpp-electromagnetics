#include "../include/sandt3D.hpp"
#include "../include/globals3d.hpp"

#include <Eigen/Dense>

#include <array>
#include <cmath>
#include <stdexcept>

namespace
{
Eigen::Vector3d cross_vec(const Eigen::Vector3d &a, const Eigen::Vector3d &b)
{
    return a.cross(b);
}
}

// Translation of sandt3D.m:
// elemental S and T matrices for 3D CT/LN (Whitney) tetrahedra.
std::pair<Eigen::MatrixXd, Eigen::MatrixXd> sandt3D(const std::vector<double> &vertex1,
                                                    const std::vector<double> &vertex2,
                                                    const std::vector<double> &vertex3,
                                                    const std::vector<double> &vertex4)
{
    const double x1 = vertex1[0], y1 = vertex1[1], z1 = vertex1[2];
    const double x2 = vertex2[0], y2 = vertex2[1], z2 = vertex2[2];
    const double x3 = vertex3[0], y3 = vertex3[1], z3 = vertex3[2];
    const double x4 = vertex4[0], y4 = vertex4[1], z4 = vertex4[2];

    Eigen::Matrix4d volume_matrix;
    volume_matrix << 1.0, x1, y1, z1,
        1.0, x2, y2, z2,
        1.0, x3, y3, z3,
        1.0, x4, y4, z4;
    const double volume = std::abs(volume_matrix.determinant()) / 6.0;
    if (volume <= 0.0)
        throw std::runtime_error("Elemento tetraedrico degenerado em sandt3D.");

    Eigen::Matrix4d A;
    A << x1, x2, x3, x4,
        y1, y2, y3, y4,
        z1, z2, z3, z4,
        1.0, 1.0, 1.0, 1.0;
    const Eigen::Matrix4d temp = A.inverse();

    std::array<Eigen::Vector3d, 4> nabla_lambda;
    Eigen::Matrix4d phi = Eigen::Matrix4d::Zero();
    std::array<std::array<Eigen::Vector3d, 4>, 4> v;

    for (int ii = 0; ii < 4; ++ii)
    {
        const double b = temp(ii, 0);
        const double c = temp(ii, 1);
        const double d = temp(ii, 2);
        nabla_lambda[ii] = {b, c, d};
    }

    for (int ii = 0; ii < 4; ++ii)
    {
        for (int jj = 0; jj < 4; ++jj)
        {
            phi(ii, jj) = nabla_lambda[ii].dot(nabla_lambda[jj]);
            v[ii][jj] = cross_vec(nabla_lambda[ii], nabla_lambda[jj]);
        }
    }

    Eigen::Matrix4d M;
    M << 2.0, 1.0, 1.0, 1.0,
        1.0, 2.0, 1.0, 1.0,
        1.0, 1.0, 2.0, 1.0,
        1.0, 1.0, 1.0, 2.0;
    M /= 20.0;

    Eigen::MatrixXd S = Eigen::MatrixXd::Zero(6, 6);
    Eigen::MatrixXd T = Eigen::MatrixXd::Zero(6, 6);
    for (int ii = 0; ii < 6; ++ii)
    {
        for (int jj = 0; jj < 6; ++jj)
        {
            const int i1 = LOCALEDGENODES[ii][0];
            const int i2 = LOCALEDGENODES[ii][1];
            const int j1 = LOCALEDGENODES[jj][0];
            const int j2 = LOCALEDGENODES[jj][1];

            S(ii, jj) = 4.0 * volume * v[i1][i2].dot(v[j1][j2]);
            T(ii, jj) = volume *
                        (phi(i2, j2) * M(i1, j1) - phi(i2, j1) * M(i1, j2) -
                         phi(i1, j2) * M(i2, j1) + phi(i1, j1) * M(i2, j2));
        }
    }

    return {S, T};
}
