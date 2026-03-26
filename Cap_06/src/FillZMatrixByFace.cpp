#include "FillZMatrixByFace.hpp"

#include "Int_pq.hpp"
#include "mom3d_globals.hpp"

namespace cap06 {

namespace {

struct FacePotentialValues {
    Eigen::Matrix<Complex, 3, 3> MagVecPot;
    Complex ScalPot;
};

FacePotentialValues Potentials(
    const int field_pt,
    const int source_pt,
    const double k,
    const Eigen::MatrixXd& r_c,
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

    FacePotentialValues values;
    for (int node = 0; node < 3; ++node) {
        values.MagVecPot.col(node) = mu_0 / (4.0 * PI) *
            (r1.cast<Complex>() * integrals.Ipq_xi + r2.cast<Complex>() * integrals.Ipq_eta +
             r3.cast<Complex>() * integrals.Ipq_zeta -
             NODE_COORD.row(qnodes(node)).transpose().cast<Complex>() * integrals.Ipq);
    }
    values.ScalPot = integrals.Ipq / (j_complex() * 2.0 * PI * omega * eps_0);
    return values;
}

}  // namespace

Eigen::MatrixXcd FillZMatrixByFace(
    const double omega,
    const double eps_0,
    const double mu_0,
    const double k,
    const Eigen::MatrixXd& r_c,
    const Eigen::MatrixXd& rho_c_pls,
    const Eigen::MatrixXd& rho_c_mns,
    const int quad_pts,
    const bool sing,
    const Eigen::VectorXi& dof2edge,
    const Eigen::VectorXi& dof_RWG)
{
    Eigen::MatrixXcd Z = Eigen::MatrixXcd::Zero(NUM_DOFS, NUM_DOFS);

    for (int melem = 0; melem < NUM_ELEMENTS; ++melem) {
        for (int nelem = 0; nelem < NUM_ELEMENTS; ++nelem) {
            const FacePotentialValues potentials =
                Potentials(melem, nelem, k, r_c, quad_pts, sing, eps_0, mu_0, omega);

            for (int nedge = 0; nedge < 3; ++nedge) {
                const int nn_id = dof_RWG(ELEMENT_EDGES(nelem, nedge));
                if (nn_id > 0) {
                    const int nn = nn_id - 1;
                    const int nnode = LOCALVERTEX(nedge);
                    const Vector3cd EllMagVecPot = ELL(dof2edge(nn)) * potentials.MagVecPot.col(nnode);
                    const Complex EllScalPot = ELL(dof2edge(nn)) * potentials.ScalPot;

                    const int source_sign = ELEMENT_PLS_MNS(nelem, nedge);
                    const Vector3cd Amn = static_cast<double>(source_sign) * EllMagVecPot;
                    const Complex Phi_mn = -static_cast<double>(source_sign) * EllScalPot;

                    for (int medge = 0; medge < 3; ++medge) {
                        const int mm_id = dof_RWG(ELEMENT_EDGES(melem, medge));
                        if (mm_id > 0) {
                            const int mm = mm_id - 1;
                            const int field_sign = ELEMENT_PLS_MNS(melem, medge);
                            const Vector3d rho_c =
                                (field_sign == +1) ? rho_c_pls.row(mm).transpose() : rho_c_mns.row(mm).transpose();

                            Z(mm, nn) += ELL(dof2edge(mm)) *
                                (j_complex() * omega * plain_dot(Amn, rho_c) / 2.0 -
                                 Phi_mn * static_cast<double>(field_sign));
                        }
                    }
                }
            }
        }
    }

    return Z;
}

}  // namespace cap06
