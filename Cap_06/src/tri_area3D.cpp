#include "tri_area3D.hpp"

namespace cap06 {

double tri_area3D(const Eigen::Vector3d& a, const Eigen::Vector3d& b, const Eigen::Vector3d& c)
{
    const Eigen::Vector3d vec_area = 0.5 * (a - c).cross(b - c);
    return vec_area.norm();
}

}  // namespace cap06
