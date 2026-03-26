#include "../include/LTQN3D.hpp"
#include "../include/globals3d.hpp"

#include <Eigen/Dense>

#include <array>

// Translation of LTQN3D.m:
// 20 LTQN basis functions on a tetrahedron.
Eigen::Matrix<double, 20, 3> LTQN3D(const std::array<double, 4> &lambda,
                                    const std::array<Eigen::Vector3d, 4> &nabla_lambda)
{
    Eigen::Matrix<double, 20, 3> functions = Eigen::Matrix<double, 20, 3>::Zero();

    for (int iedge = 0; iedge < 6; ++iedge)
    {
        const int n1 = LOCALEDGENODES[iedge][0];
        const int n2 = LOCALEDGENODES[iedge][1];

        functions.row(iedge) =
            (lambda[n1] * nabla_lambda[n2] - lambda[n2] * nabla_lambda[n1]).transpose();
        functions.row(iedge + 6) =
            (lambda[n1] * nabla_lambda[n2] + lambda[n2] * nabla_lambda[n1]).transpose();
    }

    for (int iface = 0; iface < 4; ++iface)
    {
        const int n1 = LOCALFACENODES[iface][0];
        const int n2 = LOCALFACENODES[iface][1];
        const int n3 = LOCALFACENODES[iface][2];

        functions.row(iface + 12) =
            (lambda[n2] * lambda[n3] * nabla_lambda[n1] +
             lambda[n1] * lambda[n3] * nabla_lambda[n2] -
             2.0 * lambda[n1] * lambda[n2] * nabla_lambda[n3])
                .transpose();

        functions.row(iface + 16) =
            (lambda[n3] * lambda[n1] * nabla_lambda[n2] +
             lambda[n2] * lambda[n1] * nabla_lambda[n3] -
             2.0 * lambda[n2] * lambda[n3] * nabla_lambda[n1])
                .transpose();
    }

    return functions;
}
