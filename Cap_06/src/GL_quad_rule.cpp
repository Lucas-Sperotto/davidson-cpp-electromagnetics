#include "GL_quad_rule.hpp"

#include <cmath>
#include <stdexcept>

namespace cap06 {

LineQuadrature GL_quad_rule(const int num_pts, const bool unit_length)
{
    LineQuadrature rule;

    switch (num_pts) {
    case 1:
        rule.points = Eigen::VectorXd(1);
        rule.weights = Eigen::VectorXd(1);
        rule.points << 0.0;
        rule.weights << 2.0;
        break;
    case 2:
        rule.points = Eigen::VectorXd(2);
        rule.weights = Eigen::VectorXd(2);
        rule.points << -1.0 / std::sqrt(3.0), 1.0 / std::sqrt(3.0);
        rule.weights << 1.0, 1.0;
        break;
    case 3:
        rule.points = Eigen::VectorXd(3);
        rule.weights = Eigen::VectorXd(3);
        rule.points << -std::sqrt(3.0 / 5.0), 0.0, std::sqrt(3.0 / 5.0);
        rule.weights << 5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0;
        break;
    case 4:
        rule.points = Eigen::VectorXd(4);
        rule.weights = Eigen::VectorXd(4);
        rule.points << -std::sqrt((3.0 + 2.0 * std::sqrt(6.0 / 5.0)) / 7.0),
            -std::sqrt((3.0 - 2.0 * std::sqrt(6.0 / 5.0)) / 7.0),
            std::sqrt((3.0 - 2.0 * std::sqrt(6.0 / 5.0)) / 7.0),
            std::sqrt((3.0 + 2.0 * std::sqrt(6.0 / 5.0)) / 7.0);
        rule.weights << (18.0 - std::sqrt(30.0)) / 36.0,
            (18.0 + std::sqrt(30.0)) / 36.0,
            (18.0 + std::sqrt(30.0)) / 36.0,
            (18.0 - std::sqrt(30.0)) / 36.0;
        break;
    case 5:
        rule.points = Eigen::VectorXd(5);
        rule.weights = Eigen::VectorXd(5);
        rule.points << -std::sqrt(5.0 + 2.0 * std::sqrt(10.0 / 7.0)) / 3.0,
            -std::sqrt(5.0 - 2.0 * std::sqrt(10.0 / 7.0)) / 3.0,
            0.0,
            std::sqrt(5.0 - 2.0 * std::sqrt(10.0 / 7.0)) / 3.0,
            std::sqrt(5.0 + 2.0 * std::sqrt(10.0 / 7.0)) / 3.0;
        rule.weights << (322.0 - 13.0 * std::sqrt(70.0)) / 900.0,
            (322.0 + 13.0 * std::sqrt(70.0)) / 900.0,
            128.0 / 225.0,
            (322.0 + 13.0 * std::sqrt(70.0)) / 900.0,
            (322.0 - 13.0 * std::sqrt(70.0)) / 900.0;
        break;
    default:
        throw std::runtime_error("Unimplemented quadrature rule requested in GL_quad_rule");
    }

    if (unit_length) {
        rule.points = 0.5 * rule.points.array() + 0.5;
        rule.weights *= 0.5;
    }

    return rule;
}

}  // namespace cap06
