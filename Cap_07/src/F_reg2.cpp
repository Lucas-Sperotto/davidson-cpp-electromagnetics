#include "../include/chapter7_functions.hpp"

#include <cmath>

// Translation of F_reg2.m:
// region 2 with lambda = k_0*cosh(t), subtracting the pole contribution.
std::complex<double> F_reg2(double t, const SommerfeldContext &context)
{
    const double lambda = context.k_0 * std::cosh(t);
    const std::complex<double> f_sing = context.residue / (lambda - context.lambda_p);
    return (F(lambda, context) - f_sing) * context.k_0 * std::sinh(t);
}
