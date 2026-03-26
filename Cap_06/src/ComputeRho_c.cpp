#include "ComputeRho_c.hpp"

#include "mom3d_globals.hpp"

namespace cap06 {

RhoCenters ComputeRho_c(const Eigen::MatrixXd& r_c)
{
    RhoCenters result;
    result.rho_c_pls = Eigen::MatrixXd::Zero(NUM_DOFS, 3);
    result.rho_c_mns = Eigen::MatrixXd::Zero(NUM_DOFS, 3);
    result.x_grid_pls = Eigen::VectorXd::Zero(NUM_DOFS);
    result.y_grid_pls = Eigen::VectorXd::Zero(NUM_DOFS);
    result.x_grid_mns = Eigen::VectorXd::Zero(NUM_DOFS);
    result.y_grid_mns = Eigen::VectorXd::Zero(NUM_DOFS);

    for (int mm = 0; mm < NUM_DOFS; ++mm) {
        const int pp_pls = EDGECONXELEMS(mm, 0);
        const int loc_node_pls = LOCALVERTEX(DOFLOCALNUM(mm, 0));
        const Vector3d vertx_pls = NODE_COORD.row(ELEMENTS(pp_pls, loc_node_pls)).transpose();
        result.rho_c_pls.row(mm) = r_c.row(pp_pls) - vertx_pls.transpose();
        result.x_grid_pls(mm) = r_c(pp_pls, 0);
        result.y_grid_pls(mm) = r_c(pp_pls, 1);

        const int pp_mns = EDGECONXELEMS(mm, 1);
        const int loc_node_mns = LOCALVERTEX(DOFLOCALNUM(mm, 1));
        const Vector3d vertx_mns = NODE_COORD.row(ELEMENTS(pp_mns, loc_node_mns)).transpose();
        result.rho_c_mns.row(mm) = -(r_c.row(pp_mns) - vertx_mns.transpose());
        result.x_grid_mns(mm) = r_c(pp_mns, 0);
        result.y_grid_mns(mm) = r_c(pp_mns, 1);
    }

    return result;
}

}  // namespace cap06
