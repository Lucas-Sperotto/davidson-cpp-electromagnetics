#ifndef CAP_07_CHAPTER7_FUNCTIONS_HPP
#define CAP_07_CHAPTER7_FUNCTIONS_HPP

#include "chapter7_common.hpp"

#include <complex>

std::complex<double> D_TE(const std::complex<double> &u_0,
                          const std::complex<double> &u,
                          double d);

std::complex<double> D_TM(const std::complex<double> &eps_r,
                          const std::complex<double> &u_0,
                          const std::complex<double> &u,
                          double d);

std::complex<double> F(double lambda, const SommerfeldContext &context);
std::complex<double> F_reg1(double t, const SommerfeldContext &context);
std::complex<double> F_reg2(double t, const SommerfeldContext &context);
std::complex<double> F_reg3(double lambda, const SommerfeldContext &context);
std::complex<double> F_static(double lambda, const SommerfeldContext &context);

double root_D_TM(double k_0,
                 const std::complex<double> &k,
                 const std::complex<double> &eps_r,
                 double d,
                 int max_iter);

std::complex<double> compute_residue(const SommerfeldContext &context);
VIntResult V_int_details(const SommerfeldContext &context);
std::complex<double> V_int(const SommerfeldContext &context);

#endif
