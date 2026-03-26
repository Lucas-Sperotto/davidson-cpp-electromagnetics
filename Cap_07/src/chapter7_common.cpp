#include "../include/chapter7_common.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace
{
std::complex<double> simpson_rule(const std::complex<double> &fa,
                                  const std::complex<double> &fm,
                                  const std::complex<double> &fb,
                                  double a,
                                  double b)
{
    return (b - a) * (fa + 4.0 * fm + fb) / 6.0;
}

std::complex<double> adaptive_simpson_recursive(
    const std::function<std::complex<double>(double)> &function,
    double a,
    double b,
    double tolerance,
    int depth,
    const std::complex<double> &fa,
    const std::complex<double> &fm,
    const std::complex<double> &fb,
    const std::complex<double> &whole)
{
    const double midpoint = 0.5 * (a + b);
    const double left_midpoint = 0.5 * (a + midpoint);
    const double right_midpoint = 0.5 * (midpoint + b);

    const std::complex<double> f_left_mid = function(left_midpoint);
    const std::complex<double> f_right_mid = function(right_midpoint);

    const std::complex<double> left =
        simpson_rule(fa, f_left_mid, fm, a, midpoint);
    const std::complex<double> right =
        simpson_rule(fm, f_right_mid, fb, midpoint, b);
    const std::complex<double> delta = left + right - whole;

    if (depth <= 0 || std::abs(delta) <= 15.0 * tolerance)
        return left + right + delta / 15.0;

    return adaptive_simpson_recursive(function, a, midpoint, tolerance / 2.0, depth - 1,
                                      fa, f_left_mid, fm, left) +
           adaptive_simpson_recursive(function, midpoint, b, tolerance / 2.0, depth - 1,
                                      fm, f_right_mid, fb, right);
}
} // namespace

std::complex<double> coth_complex(const std::complex<double> &z)
{
    if (std::real(z) > 20.0)
    {
        const std::complex<double> e_minus = std::exp(-z);
        return (1.0 + e_minus * e_minus) / (1.0 - e_minus * e_minus);
    }
    if (std::real(z) < -20.0)
    {
        const std::complex<double> e_plus = std::exp(z);
        return -(1.0 + e_plus * e_plus) / (1.0 - e_plus * e_plus);
    }
    return std::cosh(z) / std::sinh(z);
}

std::complex<double> sech_complex(const std::complex<double> &z)
{
    if (std::real(z) > 20.0)
    {
        const std::complex<double> e_minus = std::exp(-z);
        return 2.0 * e_minus / (1.0 + e_minus * e_minus);
    }
    if (std::real(z) < -20.0)
    {
        const std::complex<double> e_plus = std::exp(z);
        return 2.0 * e_plus / (1.0 + e_plus * e_plus);
    }
    return 1.0 / std::cosh(z);
}

std::complex<double> csch_complex(const std::complex<double> &z)
{
    if (std::real(z) > 20.0)
    {
        const std::complex<double> e_minus = std::exp(-z);
        return 2.0 * e_minus / (1.0 - e_minus * e_minus);
    }
    if (std::real(z) < -20.0)
    {
        const std::complex<double> e_plus = std::exp(z);
        return -2.0 * e_plus / (1.0 - e_plus * e_plus);
    }
    return 1.0 / std::sinh(z);
}

std::complex<double> adaptive_simpson(
    const std::function<std::complex<double>(double)> &function,
    double start,
    double stop,
    double tolerance,
    int max_depth)
{
    if (!(start <= stop))
        throw std::invalid_argument("adaptive_simpson requer intervalo crescente.");

    const double midpoint = 0.5 * (start + stop);
    const std::complex<double> fa = function(start);
    const std::complex<double> fm = function(midpoint);
    const std::complex<double> fb = function(stop);
    const std::complex<double> whole = simpson_rule(fa, fm, fb, start, stop);

    return adaptive_simpson_recursive(function, start, stop, tolerance, max_depth,
                                      fa, fm, fb, whole);
}

std::complex<double> trapz_uniform(const std::vector<std::complex<double>> &values,
                                   double delta)
{
    if (values.size() < 2)
        return 0.0;

    std::complex<double> sum = 0.5 * (values.front() + values.back());
    for (std::size_t index = 1; index + 1 < values.size(); ++index)
        sum += values[index];

    return delta * sum;
}

std::complex<double> interp1_linear(const std::vector<double> &x,
                                    const std::vector<std::complex<double>> &y,
                                    double x_query)
{
    if (x.size() != y.size() || x.empty())
        throw std::invalid_argument("interp1_linear requer vetores nao vazios de mesmo tamanho.");

    if (x_query <= x.front())
        return y.front();
    if (x_query >= x.back())
        return y.back();

    const auto upper = std::upper_bound(x.begin(), x.end(), x_query);
    const std::size_t right = static_cast<std::size_t>(upper - x.begin());
    const std::size_t left = right - 1;

    const double alpha = (x_query - x[left]) / (x[right] - x[left]);
    return (1.0 - alpha) * y[left] + alpha * y[right];
}
