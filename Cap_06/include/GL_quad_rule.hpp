#pragma once

#include <Eigen/Dense>

namespace cap06 {

struct LineQuadrature {
    Eigen::VectorXd points;
    Eigen::VectorXd weights;
};

LineQuadrature GL_quad_rule(int num_pts, bool unit_length = false);

}  // namespace cap06
