#pragma once

#include <Eigen/Dense>

namespace cap06 {

struct RhoCenters {
    Eigen::MatrixXd rho_c_pls;
    Eigen::MatrixXd rho_c_mns;
    Eigen::VectorXd x_grid_pls;
    Eigen::VectorXd y_grid_pls;
    Eigen::VectorXd x_grid_mns;
    Eigen::VectorXd y_grid_mns;
};

RhoCenters ComputeRho_c(const Eigen::MatrixXd& r_c);

}  // namespace cap06
