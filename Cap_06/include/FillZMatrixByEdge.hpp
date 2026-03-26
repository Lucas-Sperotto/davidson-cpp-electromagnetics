#pragma once

#include <Eigen/Dense>

namespace cap06 {

Eigen::MatrixXcd FillZMatrixByEdge(
    double omega,
    double eps_0,
    double mu_0,
    double k,
    const Eigen::MatrixXd& r_c,
    const Eigen::MatrixXd& rho_c_pls,
    const Eigen::MatrixXd& rho_c_mns,
    int quad_pts,
    bool sing,
    const Eigen::VectorXi& dof2edge);

}  // namespace cap06
