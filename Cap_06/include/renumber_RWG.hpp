#pragma once

#include <Eigen/Dense>

#include <utility>

namespace cap06 {

std::pair<Eigen::VectorXi, Eigen::VectorXi> renumber_RWG(const Eigen::VectorXi& dof_free_flag);

}  // namespace cap06
