// Cap_02/src/gaussder_norm.hpp
#ifndef GAUSSDER_NORM_HPP
#define GAUSSDER_NORM_HPP

#include <cmath>

inline double gaussder_norm(double t, double m, double sigma) {
    return -std::exp(0.5) * (t - m) / sigma * std::exp(-std::pow((t - m), 2) / (2 * sigma * sigma));
}

#endif // GAUSSDER_NORM_HPP
