#include "../include/curl_LTQN3D.hpp"
#include "../include/globals3d.hpp"

#include <Eigen/Dense>

#include <array>

// Translation of curl_LTQN3D.m:
// curls of the 20 LTQN basis functions.
Eigen::Matrix<double, 20, 3> curl_LTQN3D(const std::array<double, 4> &lambda,
                                         const std::array<Eigen::Vector3d, 4> &nabla_lambda)
{
    Eigen::Matrix<double, 20, 3> curls = Eigen::Matrix<double, 20, 3>::Zero();

    for (int iedge = 0; iedge < 6; ++iedge)
    {
        const int n1 = LOCALEDGENODES[iedge][0];
        const int n2 = LOCALEDGENODES[iedge][1];
        curls.row(iedge) =
            (2.0 * nabla_lambda[n1].cross(nabla_lambda[n2])).transpose();
    }

    for (int iface = 0; iface < 4; ++iface)
    {
        const int n1 = LOCALFACENODES[iface][0];
        const int n2 = LOCALFACENODES[iface][1];
        const int n3 = LOCALFACENODES[iface][2];

        curls.row(iface + 12) =
            (-3.0 * lambda[n2] * nabla_lambda[n1].cross(nabla_lambda[n3]) -
             3.0 * lambda[n1] * nabla_lambda[n2].cross(nabla_lambda[n3]))
                .transpose();

        curls.row(iface + 16) =
            (-3.0 * lambda[n3] * nabla_lambda[n2].cross(nabla_lambda[n1]) -
             3.0 * lambda[n2] * nabla_lambda[n3].cross(nabla_lambda[n1]))
                .transpose();
    }

    return curls;
}
