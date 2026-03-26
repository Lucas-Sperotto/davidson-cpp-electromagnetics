#include "intg_sing_SGF.hpp"

#include "GL_quad_rule.hpp"
#include "hat.hpp"
#include "simplex_area.hpp"

namespace cap06 {

SingularIntegral intg_sing_SGF(
    const double k,
    const Vector3d& r,
    const Vector3d& r1,
    const Vector3d& r2,
    const Vector3d& r3,
    const int num_pts_rad,
    const int num_pts_trans)
{
    Eigen::Matrix<Complex, 3, 1> I = Eigen::Matrix<Complex, 3, 1>::Zero();
    Eigen::Matrix<Complex, 3, 1> I1 = Eigen::Matrix<Complex, 3, 1>::Zero();
    Eigen::Matrix<Complex, 3, 1> I2 = Eigen::Matrix<Complex, 3, 1>::Zero();

    const Vector3d N = (r2 - r1).cross(r3 - r1);
    const Vector3d r0 = r + ((r1 - r).dot(N) / N.squaredNorm()) * N;
    const Vector3d r_vec = r - r0;
    const double z = hat(N).dot(r_vec);
    const Vector3d xi0 = simplex_area(r0, r1, r2, r3);

    for (int tt = 0; tt < 3; ++tt) {
        Vector3d r1_p = r0;
        Vector3d r2_p;
        Vector3d r3_p;
        Eigen::Matrix3d T;

        switch (tt) {
        case 0:
            r2_p = r2;
            r3_p = r3;
            T << xi0(0), 0.0, 0.0,
                xi0(1), 1.0, 0.0,
                xi0(2), 0.0, 1.0;
            break;
        case 1:
            r2_p = r3;
            r3_p = r1;
            T << xi0(0), 0.0, 1.0,
                xi0(1), 0.0, 0.0,
                xi0(2), 1.0, 0.0;
            break;
        default:
            r2_p = r1;
            r3_p = r2;
            T << xi0(0), 1.0, 0.0,
                xi0(1), 0.0, 1.0,
                xi0(2), 0.0, 0.0;
            break;
        }

        const Vector3d ell1_p = r3_p - r2_p;
        const Vector3d ell2_p = r1_p - r3_p;
        const Vector3d ell3_p = r2_p - r1_p;
        const Vector3d hat_n_p = hat(ell1_p.cross(ell2_p));
        const double A_p = 0.5 * hat_n_p.dot(ell1_p.cross(ell2_p));
        const Vector3d h1_p = (2.0 * A_p / ell1_p.squaredNorm()) * ell1_p.cross(hat_n_p);

        const LineQuadrature rule_i = GL_quad_rule(num_pts_trans, true);
        const LineQuadrature rule_j = GL_quad_rule(num_pts_rad, true);

        Eigen::Vector3d xi_pij = Eigen::Vector3d::Zero();
        for (int jj = 0; jj < num_pts_rad; ++jj) {
            xi_pij(0) = rule_j.points(jj);
            const double y_pj = h1_p.norm() * (1.0 - rule_j.points(jj));
            const double x_Lj = hat_n_p.dot(hat(h1_p).cross(ell2_p)) * (1.0 - xi_pij(0));
            const double x_Uj = -hat_n_p.dot(hat(h1_p).cross(ell3_p)) * (1.0 - xi_pij(0));
            const double root_term = std::sqrt(y_pj * y_pj + z * z);
            const double u_Lj = std::asinh(x_Lj / root_term);
            const double u_Uj = std::asinh(x_Uj / root_term);

            for (int ii = 0; ii < num_pts_trans; ++ii) {
                const double u_ij = u_Lj * (1.0 - rule_i.points(ii)) + u_Uj * rule_i.points(ii);
                const double x_pij = root_term * std::sinh(u_ij);
                xi_pij(2) =
                    hat_n_p.dot(ell3_p.cross(hat(h1_p) * y_pj - hat(ell1_p) * x_pij)) / (2.0 * A_p);
                xi_pij(1) = 1.0 - xi_pij(2) - xi_pij(0);

                const Vector3d xi_ij = T * xi_pij;
                const Vector3d mapped_r = r1 * xi_ij(0) + r2 * xi_ij(1) + r3 * xi_ij(2);
                const double R_ij = (r - mapped_r).norm();
                const Complex factor = rule_i.weights(ii) * rule_j.weights(jj) * h1_p.norm() * (u_Uj - u_Lj) *
                    std::exp(-j_complex() * k * R_ij);

                I(tt) += factor;
                I1(tt) += factor * xi_ij(0);
                I2(tt) += factor * xi_ij(1);
            }
        }
    }

    SingularIntegral result;
    result.value = I.sum();
    result.xi1 = I1.sum();
    result.xi2 = I2.sum();
    result.xi3 = result.value - result.xi1 - result.xi2;
    return result;
}

}  // namespace cap06
