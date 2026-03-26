#include "FillVVector.hpp"

#include "mom3d_globals.hpp"

#include <stdexcept>

namespace cap06 {

Eigen::VectorXcd FillVVector(
    const Eigen::MatrixXd& rho_c_pls,
    const Eigen::MatrixXd& rho_c_mns,
    const double EMag,
    const double theta_0,
    const double phi_0,
    const Eigen::VectorXi& dof2edge)
{
    if (theta_0 != 0.0 || phi_0 != 0.0) {
        throw std::runtime_error("Only normally incident plane wave implemented at present");
    }

    Eigen::VectorXcd V = Eigen::VectorXcd::Zero(NUM_DOFS);
    Eigen::MatrixXd E = Eigen::MatrixXd::Zero(NUM_ELEMENTS, 3);
    E.col(0).setConstant(EMag);

    for (int mm = 0; mm < NUM_DOFS; ++mm) {
        const int pp_pls = EDGECONXELEMS(mm, 0);
        const int pp_mns = EDGECONXELEMS(mm, 1);
        const double sampled = 0.5 * E.row(pp_pls).dot(rho_c_pls.row(mm)) +
            0.5 * E.row(pp_mns).dot(rho_c_mns.row(mm));
        V(mm) = ELL(dof2edge(mm)) * sampled;
    }

    return V;
}

}  // namespace cap06
