#include "../include/chapter7_functions.hpp"

#include <cmath>
#include <numbers>
#include <stdexcept>
#include <string>

namespace
{
bool is_finite_complex(const std::complex<double> &value)
{
    return std::isfinite(value.real()) && std::isfinite(value.imag());
}
} // namespace

std::complex<double> compute_residue(const SommerfeldContext &context)
{
    const std::complex<double> u =
        std::sqrt(std::complex<double>(context.lambda_p * context.lambda_p, 0.0) - context.k * context.k);
    const std::complex<double> u_0 =
        std::sqrt(std::complex<double>(context.lambda_p * context.lambda_p - context.k_0 * context.k_0, 0.0));

    const std::complex<double> dD_TE_by_dlambda =
        context.lambda_p / u_0 +
        (context.lambda_p / u) * coth_complex(u * context.d) -
        context.lambda_p * context.d * std::pow(csch_complex(u * context.d), 2);

    const std::complex<double> dD_TM_by_dlambda =
        context.eps_r * context.lambda_p / u_0 +
        (context.lambda_p / u) * std::tanh(u * context.d) +
        context.lambda_p * context.d * std::pow(sech_complex(u * context.d), 2);

    const std::complex<double> numerator =
        std::cyl_bessel_j(0.0, context.lambda_p * context.rho) * context.lambda_p *
        (u_0 + u * std::tanh(u * context.d));

    const std::complex<double> denominator =
        D_TM(context.eps_r, u_0, u, context.d) * dD_TE_by_dlambda +
        D_TE(u_0, u, context.d) * dD_TM_by_dlambda;

    const std::complex<double> residue = numerator / denominator;
    if (!is_finite_complex(residue))
    {
        throw std::runtime_error(
            "compute_residue produziu NaN/Inf para rho = " + std::to_string(context.rho));
    }
    return residue;
}

VIntResult V_int_details(const SommerfeldContext &context)
{
    if (context.lambda_p <= 0.0)
        throw std::invalid_argument("V_int_details requer lambda_p positivo.");

    VIntResult result;
    result.lambda_p = context.lambda_p;
    result.residue = compute_residue(context);

    SommerfeldContext work = context;
    work.residue = result.residue;

    result.reg1 = adaptive_simpson(
        [&](double t)
        { return F_reg1(t, work); },
        0.0, std::numbers::pi / 2.0);
    if (!is_finite_complex(result.reg1))
        throw std::runtime_error("V_int: regiao 1 produziu NaN/Inf.");

    const double zeta = 0.999;
    const double region2_stop = std::acosh(zeta * std::sqrt(context.eps_r_prime));
    result.reg2_xsing = adaptive_simpson(
        [&](double t)
        { return F_reg2(t, work); },
        0.0, region2_stop);
    if (!is_finite_complex(result.reg2_xsing))
        throw std::runtime_error("V_int: regiao 2 sem singularidade produziu NaN/Inf.");

    result.reg2_sing =
        result.residue *
            std::log((context.k_0 * std::sqrt(context.eps_r_prime) - context.lambda_p) /
                     (context.lambda_p - context.k_0)) -
        std::complex<double>(0.0, std::numbers::pi) * result.residue;
    if (!is_finite_complex(result.reg2_sing))
        throw std::runtime_error("V_int: termo singular da regiao 2 produziu NaN/Inf.");

    result.reg2 = result.reg2_xsing + result.reg2_sing;
    if (!is_finite_complex(result.reg2))
        throw std::runtime_error("V_int: regiao 2 produziu NaN/Inf.");

    const double region3_start =
        (1.0 / zeta) * context.k_0 * std::sqrt(context.eps_r_prime);
    double lambda_upper = 10.0 * context.k_0;

    std::complex<double> region3_xstatic1 = adaptive_simpson(
        [&](double lambda)
        { return F_reg3(lambda, work); },
        region3_start, lambda_upper);
    if (!is_finite_complex(region3_xstatic1))
        throw std::runtime_error("V_int: primeira integral da regiao 3 produziu NaN/Inf.");
    std::complex<double> region3_xstatic2 = adaptive_simpson(
        [&](double lambda)
        { return F_reg3(lambda, work); },
        region3_start, 2.0 * lambda_upper);
    if (!is_finite_complex(region3_xstatic2))
        throw std::runtime_error("V_int: segunda integral da regiao 3 produziu NaN/Inf.");

    std::complex<double> err = 0.0;
    if (std::abs(region3_xstatic2) > 0.0)
        err = (region3_xstatic1 - region3_xstatic2) / std::abs(region3_xstatic2);
    else
        err = region3_xstatic1 - region3_xstatic2;

    const double eps_rel = std::max(std::abs(result.reg2) * 1e-4, 1e-8);
    int count = 0;
    const int max_count = 10;
    while (std::abs(err) > eps_rel && count < max_count)
    {
        ++count;
        region3_xstatic1 = region3_xstatic2;
        lambda_upper *= 2.0;
        region3_xstatic2 = adaptive_simpson(
            [&](double lambda)
            { return F_reg3(lambda, work); },
            region3_start, 2.0 * lambda_upper);
        if (!is_finite_complex(region3_xstatic2))
            throw std::runtime_error("V_int: refinamento da regiao 3 produziu NaN/Inf.");

        if (std::abs(region3_xstatic2) > 0.0)
            err = (region3_xstatic1 - region3_xstatic2) / std::abs(region3_xstatic2);
        else
            err = region3_xstatic1 - region3_xstatic2;
    }

    result.reg3_iterations = count;
    result.lambda_upper_final = 2.0 * lambda_upper;
    result.reg3_xstatic = region3_xstatic2;

    const std::complex<double> static_integral = adaptive_simpson(
        [&](double lambda)
        { return F_static(lambda, work); },
        0.0, std::sqrt(context.eps_r_prime) * context.k_0);
    if (!is_finite_complex(static_integral))
        throw std::runtime_error("V_int: integral estatica produziu NaN/Inf.");

    result.reg3 =
        result.reg3_xstatic + 1.0 / (context.rho * (1.0 + context.eps_r)) - static_integral;
    if (!is_finite_complex(result.reg3))
        throw std::runtime_error("V_int: regiao 3 produziu NaN/Inf.");

    result.value = result.reg1 + result.reg2 + result.reg3;
    if (!is_finite_complex(result.value))
        throw std::runtime_error("V_int: resultado final produziu NaN/Inf.");
    return result;
}

std::complex<double> V_int(const SommerfeldContext &context)
{
    return V_int_details(context).value;
}
