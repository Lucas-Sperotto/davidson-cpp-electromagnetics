#include "FillZMatrixByEdge.hpp"

#include "Int_pq.hpp"
#include "mom3d_globals.hpp"

namespace cap06 {

namespace {

struct PotentialValues {
    Vector3cd MagVecPot;
    Complex ScalPot;
};

PotentialValues Potentials(
    const int field_pt,
    const int source_pt,
    const int source_edge,
    const int source_tri,
    const double k,
    const Eigen::MatrixXd& r_c,
    const Eigen::VectorXi& dof2edge,
    const int quad_pts,
    const bool sing,
    const double eps_0,
    const double mu_0,
    const double omega)
{
    const IntegralPQ integrals = Int_pq(field_pt, source_pt, r_c.row(field_pt).transpose(), k, quad_pts, sing);
    const Eigen::Vector3i qnodes = ELEMENTS.row(source_pt).transpose();
    const Vector3d r1 = NODE_COORD.row(qnodes(0)).transpose();
    const Vector3d r2 = NODE_COORD.row(qnodes(1)).transpose();
    const Vector3d r3 = NODE_COORD.row(qnodes(2)).transpose();
    const int ii = DOFLOCALNUM(source_edge, source_tri);

    PotentialValues values;
    values.MagVecPot = mu_0 * ELL(dof2edge(source_edge)) / (4.0 * PI) *
        (r1.cast<Complex>() * integrals.Ipq_xi + r2.cast<Complex>() * integrals.Ipq_eta +
         r3.cast<Complex>() * integrals.Ipq_zeta -
         NODE_COORD.row(qnodes(ii)).transpose().cast<Complex>() * integrals.Ipq);
    values.ScalPot = ELL(dof2edge(source_edge)) / (j_complex() * 2.0 * PI * omega * eps_0) * integrals.Ipq;
    return values;
}

}  // namespace

Eigen::MatrixXcd FillZMatrixByEdge(
    const double omega,
    const double eps_0,
    const double mu_0,
    const double k,
    const Eigen::MatrixXd& r_c,
    const Eigen::MatrixXd& rho_c_pls,
    const Eigen::MatrixXd& rho_c_mns,
    const int quad_pts,
    const bool sing,
    const Eigen::VectorXi& dof2edge)
{
    Eigen::MatrixXcd Z = Eigen::MatrixXcd::Zero(NUM_DOFS, NUM_DOFS);

    for (int mm = 0; mm < NUM_DOFS; ++mm) {
        for (int nn = 0; nn < NUM_DOFS; ++nn) {
            const int pp_pls = EDGECONXELEMS(mm, 0);
            const int pp_mns = EDGECONXELEMS(mm, 1);
            const int qq_pls = EDGECONXELEMS(nn, 0);
            const int qq_mns = EDGECONXELEMS(nn, 1);

            const PotentialValues p_pls_pls =
                Potentials(pp_pls, qq_pls, nn, 0, k, r_c, dof2edge, quad_pts, sing, eps_0, mu_0, omega);
            const Vector3cd Amn_pls_source_pls = p_pls_pls.MagVecPot;
            const Complex Phi_mn_pls_source_pls = -p_pls_pls.ScalPot;

            const PotentialValues p_pls_mns =
                Potentials(pp_pls, qq_mns, nn, 1, k, r_c, dof2edge, quad_pts, sing, eps_0, mu_0, omega);
            const Vector3cd Amn_pls_source_mns = -p_pls_mns.MagVecPot;
            const Complex Phi_mn_pls_source_mns = p_pls_mns.ScalPot;

            const Vector3cd Amn_pls = Amn_pls_source_pls + Amn_pls_source_mns;
            const Complex Phi_mn_pls = Phi_mn_pls_source_pls + Phi_mn_pls_source_mns;

            const PotentialValues p_mns_pls =
                Potentials(pp_mns, qq_pls, nn, 0, k, r_c, dof2edge, quad_pts, sing, eps_0, mu_0, omega);
            const Vector3cd Amn_mns_source_pls = p_mns_pls.MagVecPot;
            const Complex Phi_mn_mns_source_pls = -p_mns_pls.ScalPot;

            const PotentialValues p_mns_mns =
                Potentials(pp_mns, qq_mns, nn, 1, k, r_c, dof2edge, quad_pts, sing, eps_0, mu_0, omega);
            const Vector3cd Amn_mns_source_mns = -p_mns_mns.MagVecPot;
            const Complex Phi_mn_mns_source_mns = p_mns_mns.ScalPot;

            const Vector3cd Amn_mns = Amn_mns_source_pls + Amn_mns_source_mns;
            const Complex Phi_mn_mns = Phi_mn_mns_source_pls + Phi_mn_mns_source_mns;

            const Vector3d rho_pls = rho_c_pls.row(mm).transpose();
            const Vector3d rho_mns = rho_c_mns.row(mm).transpose();

            Z(mm, nn) = j_complex() * omega * (plain_dot(Amn_pls, rho_pls) / 2.0 + plain_dot(Amn_mns, rho_mns) / 2.0) +
                Phi_mn_mns - Phi_mn_pls;
            Z(mm, nn) *= ELL(dof2edge(mm));
        }
    }

    return Z;
}

}  // namespace cap06
