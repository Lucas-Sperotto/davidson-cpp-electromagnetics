#include "../include/chapter7_functions.hpp"

#include <cmath>

// Translation of F_static.m:
// extracted static contribution in region 3.
std::complex<double> F_static(double lambda, const SommerfeldContext &context)
{
    return std::cyl_bessel_j(0.0, lambda * context.rho) / (1.0 + context.eps_r);
}
