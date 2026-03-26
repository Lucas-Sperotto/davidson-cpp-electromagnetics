#include "../include/chapter7_functions.hpp"

// Translation of D_TM.m:
// TM surface-wave denominator for the grounded substrate.
std::complex<double> D_TM(const std::complex<double> &eps_r,
                          const std::complex<double> &u_0,
                          const std::complex<double> &u,
                          double d)
{
    return eps_r * u_0 + u * std::tanh(u * d);
}
