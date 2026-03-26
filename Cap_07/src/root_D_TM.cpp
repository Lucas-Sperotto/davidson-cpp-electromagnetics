#include "../include/chapter7_functions.hpp"

#include <cmath>

// Translation of root_D_TM.m:
// locate the TM pole using bisection on the real axis.
double root_D_TM(double k_0,
                 const std::complex<double> &k,
                 const std::complex<double> &eps_r,
                 double d,
                 int max_iter)
{
    double x2 = k_0 * std::sqrt(std::real(eps_r));
    std::complex<double> u = std::sqrt(std::complex<double>(x2 * x2, 0.0) - k * k);
    std::complex<double> u_0 =
        std::sqrt(std::complex<double>(x2 * x2 - k_0 * k_0, 0.0));
    double fmid = std::real(D_TM(eps_r, u_0, u, d));

    double x1 = k_0;
    u = std::sqrt(std::complex<double>(x1 * x1, 0.0) - k * k);
    u_0 = std::sqrt(std::complex<double>(x1 * x1 - k_0 * k_0, 0.0));
    const double f = std::real(D_TM(eps_r, u_0, u, d));

    double tm_root = 0.0;
    double delta_x = 0.0;
    if (f < 0.0)
    {
        tm_root = x1;
        delta_x = x2 - x1;
    }
    else
    {
        tm_root = x2;
        delta_x = x1 - x2;
    }

    for (int kk = 0; kk < max_iter; ++kk)
    {
        delta_x /= 2.0;
        const double xmid = tm_root + delta_x;
        u = std::sqrt(std::complex<double>(xmid * xmid, 0.0) - k * k);
        u_0 = std::sqrt(std::complex<double>(xmid * xmid - k_0 * k_0, 0.0));
        fmid = std::real(D_TM(eps_r, u_0, u, d));
        if (fmid < 0.0)
            tm_root = xmid;
    }

    return tm_root;
}
