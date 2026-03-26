#pragma once

#include <Eigen/Dense>

namespace cap06 {

Eigen::VectorXcd FillVVector(
    const Eigen::MatrixXd& rho_c_pls,
    const Eigen::MatrixXd& rho_c_mns,
    double EMag,
    double theta_0,
    double phi_0,
    const Eigen::VectorXi& dof2edge);

}  // namespace cap06
