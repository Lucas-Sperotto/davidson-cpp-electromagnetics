#ifndef CAP_07_CHAPTER7_COMMON_HPP
#define CAP_07_CHAPTER7_COMMON_HPP

#include <complex>
#include <functional>
#include <vector>

struct SommerfeldContext
{
    double k_0 = 0.0;
    std::complex<double> k = 0.0;
    std::complex<double> eps_r = 1.0;
    double eps_r_prime = 1.0;
    double d = 0.0;
    double rho = 0.0;
    double lambda_p = 0.0;
    std::complex<double> residue = 0.0;
};

struct VIntResult
{
    std::complex<double> value = 0.0;
    std::complex<double> reg1 = 0.0;
    std::complex<double> reg2 = 0.0;
    std::complex<double> reg2_xsing = 0.0;
    std::complex<double> reg2_sing = 0.0;
    std::complex<double> reg3 = 0.0;
    std::complex<double> reg3_xstatic = 0.0;
    std::complex<double> residue = 0.0;
    double lambda_p = 0.0;
    double lambda_upper_final = 0.0;
    int reg3_iterations = 0;
};

std::complex<double> coth_complex(const std::complex<double> &z);
std::complex<double> sech_complex(const std::complex<double> &z);
std::complex<double> csch_complex(const std::complex<double> &z);

std::complex<double> adaptive_simpson(
    const std::function<std::complex<double>(double)> &function,
    double start,
    double stop,
    double tolerance = 1e-6,
    int max_depth = 18);

std::complex<double> trapz_uniform(const std::vector<std::complex<double>> &values,
                                   double delta = 1.0);

std::complex<double> interp1_linear(const std::vector<double> &x,
                                    const std::vector<std::complex<double>> &y,
                                    double x_query);

#endif
