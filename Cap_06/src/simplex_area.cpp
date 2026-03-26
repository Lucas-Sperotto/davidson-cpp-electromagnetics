#include "simplex_area.hpp"

#include <cmath>

namespace cap06 {

Eigen::Vector3d simplex_area(
    const Eigen::Vector3d& q,
    const Eigen::Vector3d& a,
    const Eigen::Vector3d& b,
    const Eigen::Vector3d& c)
{
    const Eigen::Vector3d ap = a - c;
    const Eigen::Vector3d bp = b - c;
    const Eigen::Vector3d qp = q - c;

    const double denom = ap.dot(ap) * bp.dot(bp) - std::pow(ap.dot(bp), 2);
    const double alpha = (bp.dot(bp) * ap.dot(qp) - ap.dot(bp) * bp.dot(qp)) / denom;
    const double beta = (ap.dot(ap) * bp.dot(qp) - ap.dot(bp) * ap.dot(qp)) / denom;

    return Eigen::Vector3d(alpha, beta, 1.0 - alpha - beta);
}

}  // namespace cap06
