#include "../include/sandt3D_LTQN.hpp"
#include "../include/LTQN3D.hpp"
#include "../include/curl_LTQN3D.hpp"
#include "../include/tet_quad.hpp"

#include <Eigen/Dense>

#include <array>
#include <cmath>
#include <stdexcept>

// Translation of sandt3D_LTQN.m:
// elemental S and T matrices for 3D LT/QN tetrahedra.
std::pair<Eigen::MatrixXd, Eigen::MatrixXd> sandt3D_LTQN(const std::vector<double> &vertex1,
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
        throw std::runtime_error("Elemento tetraedrico degenerado em sandt3D_LTQN.");

    Eigen::Matrix4d A;
    A << x1, x2, x3, x4,
        y1, y2, y3, y4,
        z1, z2, z3, z4,
        1.0, 1.0, 1.0, 1.0;
    const Eigen::Matrix4d temp = A.inverse();

    std::array<Eigen::Vector3d, 4> nabla_lambda;
    for (int ii = 0; ii < 4; ++ii)
        nabla_lambda[ii] = {temp(ii, 0), temp(ii, 1), temp(ii, 2)};

    Eigen::MatrixXd S = Eigen::MatrixXd::Zero(20, 20);
    Eigen::MatrixXd T = Eigen::MatrixXd::Zero(20, 20);

    const auto [w, lambda] = tet_quad(11);
    for (int nn = 0; nn < w.size(); ++nn)
    {
        std::array<double, 4> bary = {
            lambda(nn, 0), lambda(nn, 1), lambda(nn, 2), lambda(nn, 3)};
        const Eigen::Matrix<double, 20, 3> funcs = LTQN3D(bary, nabla_lambda);
        const Eigen::Matrix<double, 20, 3> curls = curl_LTQN3D(bary, nabla_lambda);

        for (int ii = 0; ii < 20; ++ii)
        {
            for (int jj = 0; jj < 20; ++jj)
            {
                S(ii, jj) += w(nn) * curls.row(ii).dot(curls.row(jj));
                T(ii, jj) += w(nn) * funcs.row(ii).dot(funcs.row(jj));
            }
        }
    }

    S *= volume;
    T *= volume;
    return {S, T};
}
