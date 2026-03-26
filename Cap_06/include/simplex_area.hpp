#pragma once

#include <Eigen/Dense>

namespace cap06 {

Eigen::Vector3d simplex_area(
    const Eigen::Vector3d& q,
    const Eigen::Vector3d& a,
    const Eigen::Vector3d& b,
    const Eigen::Vector3d& c);

}  // namespace cap06
