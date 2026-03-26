#ifndef CAP_11_LTQN3D_HPP
#define CAP_11_LTQN3D_HPP

#include <Eigen/Dense>

#include <array>

Eigen::Matrix<double, 20, 3> LTQN3D(const std::array<double, 4> &lambda,
                                    const std::array<Eigen::Vector3d, 4> &nabla_lambda);

#endif
