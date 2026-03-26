#include "sphereRCS.hpp"

#include "mom3d_globals.hpp"

#include <cmath>
#include <complex>

namespace cap06 {

namespace {

std::complex<double> riccati_hankel2(const int n, const double x)
{
    const double nu = static_cast<double>(n) + 0.5;
    const std::complex<double> H_nu(
        std::cyl_bessel_j(nu, x),
        -std::cyl_neumann(nu, x));
    return std::sqrt(PI * x / 2.0) * H_nu;
}

std::complex<double> riccati_hankel2_prime(const int n, const double x)
{
    return riccati_hankel2(n - 1, x) - (static_cast<double>(n) / x) * riccati_hankel2(n, x);
}

double sphereRCS_single(const double a, const double k, const int N)
{
    const double lambda = 2.0 * PI / k;
    const double ka = k * a;
    std::complex<double> rcs_sum(0.0, 0.0);

    for (int n = 1; n <= N; ++n) {
        const std::complex<double> Hcaret2n = riccati_hankel2(n, ka);
        const std::complex<double> Hcaret2nprime = riccati_hankel2_prime(n, ka);
        const double sign = (n % 2 == 0) ? 1.0 : -1.0;
        rcs_sum += sign * (2.0 * n + 1.0) / (Hcaret2n * Hcaret2nprime);
    }

    return (lambda * lambda) / (4.0 * PI) * std::norm(rcs_sum);
}

}  // namespace

std::vector<double> sphereRCS(const double a, const std::vector<double>& k, const int N)
{
    std::vector<double> result(k.size(), 0.0);
    for (std::size_t kk = 0; kk < k.size(); ++kk) {
        result[kk] = sphereRCS_single(a, k[kk], N);
    }
    return result;
}

}  // namespace cap06
