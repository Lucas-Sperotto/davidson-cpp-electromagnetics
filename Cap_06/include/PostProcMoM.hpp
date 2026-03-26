#pragma once

#include <Eigen/Dense>

namespace cap06 {

void PostProcMoM(
    const Eigen::VectorXcd& I,
    double EMag,
    const Eigen::VectorXi& dof2edge,
    double eta_0,
    double L,
    double W,
    int Xmesh,
    int Ymesh,
    int ProbType,
    int quad_pts,
    bool sing);

}  // namespace cap06
