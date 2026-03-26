// src/gaussder.cpp
#include "../include/gaussder.hpp"

#include <cmath>

double gaussder_norm(double t, double m, double sigma)
{
    return -(std::exp(0.5) / sigma) * (t - m) *
           std::exp(-std::pow(t - m, 2.0) / (2.0 * sigma * sigma));
}
