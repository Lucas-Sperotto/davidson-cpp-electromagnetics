#include "../include/chapter7_functions.hpp"

// Translation of D_TE.m:
// TE surface-wave denominator for the grounded substrate.
std::complex<double> D_TE(const std::complex<double> &u_0,
                          const std::complex<double> &u,
                          double d)
{
    return u_0 + u * coth_complex(u * d);
}
