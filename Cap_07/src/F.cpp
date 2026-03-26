#include "../include/chapter7_functions.hpp"

#include <cmath>

// Translation of F.m:
// integrand of the scalar potential V, without the 1/(2*pi*eps_0) factor.
std::complex<double> F(double lambda, const SommerfeldContext &context)
{
    const std::complex<double> u =
        std::sqrt(std::complex<double>(lambda * lambda, 0.0) - context.k * context.k);
    const std::complex<double> u_0 =
        std::sqrt(std::complex<double>(lambda * lambda - context.k_0 * context.k_0, 0.0));

    return std::cyl_bessel_j(0.0, lambda * context.rho) * lambda *
           (u_0 + u * std::tanh(u * context.d)) /
           (D_TE(u_0, u, context.d) * D_TM(context.eps_r, u_0, u, context.d));
}
