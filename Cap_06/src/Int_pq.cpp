#include "Int_pq.hpp"

#include "intg_sing_SGF.hpp"
#include "tri_area3D.hpp"
#include "tri_quad.hpp"

namespace cap06 {

IntegralPQ Int_pq(const int p, const int q, const Vector3d& r_cp, const double k, const int quad_pts, const bool sing)
{
    const Eigen::Vector3i qnodes = ELEMENTS.row(q).transpose();
    const Vector3d n1 = NODE_COORD.row(qnodes(0)).transpose();
    const Vector3d n2 = NODE_COORD.row(qnodes(1)).transpose();
    const Vector3d n3 = NODE_COORD.row(qnodes(2)).transpose();
    const double area = tri_area3D(n1, n2, n3);

    IntegralPQ result{Complex(0.0, 0.0), Complex(0.0, 0.0), Complex(0.0, 0.0), Complex(0.0, 0.0)};

    if (p == q && sing) {
        const SingularIntegral singular = intg_sing_SGF(k, r_cp, n1, n2, n3, 3, 4);
        result.Ipq = singular.value / (2.0 * area);
        result.Ipq_xi = singular.xi1 / (2.0 * area);
        result.Ipq_eta = singular.xi2 / (2.0 * area);
        result.Ipq_zeta = singular.xi3 / (2.0 * area);
        return result;
    }

    const TriangleQuadrature quad = tri_quad(quad_pts);
    const Eigen::VectorXd w = quad.weights / 2.0;

    for (int nn = 0; nn < quad_pts; ++nn) {
        const Vector3d r_prime =
            quad.lambda(nn, 0) * n1 + quad.lambda(nn, 1) * n2 + quad.lambda(nn, 2) * n3;
        const double R_p = (r_cp - r_prime).norm();
        const Complex GF = std::exp(-j_complex() * k * R_p) / R_p;

        result.Ipq += w(nn) * GF;
        result.Ipq_xi += w(nn) * quad.lambda(nn, 0) * GF;
        result.Ipq_eta += w(nn) * quad.lambda(nn, 1) * GF;
    }

    result.Ipq_zeta = result.Ipq - result.Ipq_xi - result.Ipq_eta;
    return result;
}

}  // namespace cap06
