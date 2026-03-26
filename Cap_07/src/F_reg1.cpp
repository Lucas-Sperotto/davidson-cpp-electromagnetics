#include "../include/chapter7_functions.hpp"

#include <cmath>

// Translation of F_reg1.m:
// region 1 with lambda = k_0*cos(t).
std::complex<double> F_reg1(double t, const SommerfeldContext &context)
{
    const double lambda = context.k_0 * std::cos(t);
    return F(lambda, context) * context.k_0 * std::sin(t);
}
