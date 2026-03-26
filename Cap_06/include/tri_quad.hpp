#pragma once

#include <Eigen/Dense>

namespace cap06 {

struct TriangleQuadrature {
    Eigen::VectorXd weights;
    Eigen::MatrixXd lambda;
};

TriangleQuadrature tri_quad(int n);

}  // namespace cap06
