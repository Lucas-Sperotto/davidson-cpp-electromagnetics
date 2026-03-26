#ifndef CAP_11_SANDT3D_LTQN_HPP
#define CAP_11_SANDT3D_LTQN_HPP

#include <Eigen/Dense>

#include <utility>
#include <vector>

std::pair<Eigen::MatrixXd, Eigen::MatrixXd> sandt3D_LTQN(const std::vector<double> &vertex1,
                                                         const std::vector<double> &vertex2,
                                                         const std::vector<double> &vertex3,
                                                         const std::vector<double> &vertex4);

#endif
