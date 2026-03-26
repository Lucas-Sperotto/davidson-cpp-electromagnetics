#include "hat.hpp"

#include <limits>

namespace cap06 {

Eigen::Vector3d hat(const Eigen::Vector3d& a)
{
    if (a.cwiseAbs().maxCoeff() < std::numeric_limits<double>::epsilon()) {
        return Eigen::Vector3d::Zero();
    }

    return a / a.norm();
}

}  // namespace cap06
