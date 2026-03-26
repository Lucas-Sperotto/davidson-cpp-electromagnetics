#include "../include/chapter7_functions.hpp"

// Translation of F_reg3.m:
// region 3 with the static term extracted.
std::complex<double> F_reg3(double lambda, const SommerfeldContext &context)
{
    return F(lambda, context) - F_static(lambda, context);
}
