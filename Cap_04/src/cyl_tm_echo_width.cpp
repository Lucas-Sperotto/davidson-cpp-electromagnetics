#include "../include/cyl_tm_echo_width.hpp"

#include <cmath>
#include <complex>
#include <numbers>
#include <vector>

namespace
{
std::complex<double> hankel2_n(int n, double x)
{
    return {std::cyl_bessel_j(static_cast<double>(n), x),
            -std::cyl_neumann(static_cast<double>(n), x)};
}
} // namespace

std::vector<double> cyl_tm_echo_width(double a,
                                      const std::vector<double> &k,
                                      int N,
                                      double phi)
{
    // Translation of cyl_TM_echo_width.m:
    // analytical TM echo width of a PEC circular cylinder.
    std::vector<double> lambda(k.size(), 0.0);
    std::vector<std::complex<double>> echo(k.size(), 0.0);
    for (std::size_t kk = 0; kk < k.size(); ++kk)
        lambda[kk] = 2.0 * std::numbers::pi / k[kk];

    for (int n = 0; n <= N; ++n)
    {
        const double epsilon_n = (n == 0) ? 1.0 : 2.0;
        for (std::size_t kk = 0; kk < k.size(); ++kk)
        {
            const double ka = k[kk] * a;
            const std::complex<double> ratio =
                std::cyl_bessel_j(static_cast<double>(n), ka) / hankel2_n(n, ka);
            echo[kk] += epsilon_n * ratio * std::cos(n * phi);
        }
    }

    std::vector<double> result(k.size(), 0.0);
    for (std::size_t kk = 0; kk < k.size(); ++kk)
        result[kk] = 2.0 * lambda[kk] / std::numbers::pi * std::norm(echo[kk]);
    return result;
}
