#include "include/chapter7_functions.hpp"

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numbers>
#include <vector>

namespace fs = std::filesystem;

int main()
{
    // Translation of scalar_pot.m:
    // plot the scalar-potential integrand and its three special regions.
    const fs::path out_dir = PROJECT_OUT_DIR;
    fs::create_directories(out_dir);

    const double c = 2.997925e8;
    const double freq = 10e9;
    const double lambda_0 = c / freq;
    const double eps_r_prime = 5.0;
    const double tan_delta = 0.0;
    const std::complex<double> eps_r =
        eps_r_prime * (1.0 - std::complex<double>(0.0, tan_delta));
    const double k_0 = 2.0 * std::numbers::pi * freq / c;
    const std::complex<double> k = k_0 * std::sqrt(eps_r);
    const double d = 0.2 * std::numbers::pi / k_0;
    const double rho = 3.0 / k_0;

    SommerfeldContext context;
    context.k_0 = k_0;
    context.k = k;
    context.eps_r = eps_r;
    context.eps_r_prime = eps_r_prime;
    context.d = d;
    context.rho = rho;
    context.lambda_p = root_D_TM(k_0, k, eps_r, d, 50);
    context.residue = compute_residue(context);

    std::ofstream dtm(out_dir / "scalar_pot_dtm.csv");
    dtm << "lambda_over_k0,re_d_tm,im_d_tm\n";
    dtm << std::scientific << std::setprecision(10);

    std::ofstream dte(out_dir / "scalar_pot_dte.csv");
    dte << "lambda_over_k0,abs_d_te\n";
    dte << std::scientific << std::setprecision(10);

    std::ofstream integrand(out_dir / "scalar_pot_integrand.csv");
    integrand << "lambda_over_k0,re_integrand,im_integrand\n";
    integrand << std::scientific << std::setprecision(10);

    for (int index = 0; index <= 10000; ++index)
    {
        const double lambda = (static_cast<double>(index) * 0.0005) * k_0;
        const std::complex<double> u =
            std::sqrt(std::complex<double>(lambda * lambda, 0.0) - k * k);
        const std::complex<double> u_0 =
            std::sqrt(std::complex<double>(lambda * lambda - k_0 * k_0, 0.0));
        const std::complex<double> d_tm = D_TM(eps_r, u_0, u, d);
        const std::complex<double> d_te = D_TE(u_0, u, d);
        const std::complex<double> value = F(lambda, context);

        dtm << lambda / k_0 << "," << d_tm.real() << "," << d_tm.imag() << "\n";
        dte << lambda / k_0 << "," << std::abs(d_te) << "\n";
        integrand << lambda / k_0 << "," << value.real() << "," << value.imag() << "\n";
    }

    std::ofstream region1(out_dir / "scalar_pot_region1.csv");
    region1 << "t,re_integrand,im_integrand\n";
    region1 << std::scientific << std::setprecision(10);
    for (int index = 0; index <= 157; ++index)
    {
        const double t = static_cast<double>(index) * 0.01;
        const std::complex<double> value = F_reg1(t, context);
        region1 << t << "," << value.real() << "," << value.imag() << "\n";
    }

    std::ofstream residue_scan(out_dir / "scalar_pot_residue_scan.csv");
    residue_scan << "lambda_over_lambda_p,approx_residue\n";
    residue_scan << std::scientific << std::setprecision(10);
    for (int index = 0; index <= 20; ++index)
    {
        const double lambda =
            (0.9 + 0.01 * static_cast<double>(index)) * context.lambda_p;
        const std::complex<double> value = F(lambda, context) * (lambda - context.lambda_p);
        residue_scan << lambda / context.lambda_p << "," << value.real() << "\n";
    }

    std::ofstream region2_lambda(out_dir / "scalar_pot_region2_lambda.csv");
    region2_lambda << "lambda_over_k0,re_full,re_subtracted\n";
    region2_lambda << std::scientific << std::setprecision(10);
    for (int index = 0; index <= static_cast<int>((std::sqrt(eps_r_prime) - 1.0) / 0.01); ++index)
    {
        const double lambda =
            k_0 * (1.0 + 0.01 * static_cast<double>(index));
        const std::complex<double> f_full = F(lambda, context);
        const std::complex<double> f_sing =
            context.residue / (lambda - context.lambda_p);
        region2_lambda << lambda / k_0 << "," << f_full.real() << ","
                       << (f_full - f_sing).real() << "\n";
    }

    std::ofstream region2_t(out_dir / "scalar_pot_region2_t.csv");
    region2_t << "t,abs_f,abs_f_sing,re_integrand_normalized\n";
    region2_t << std::scientific << std::setprecision(10);

    std::vector<std::complex<double>> reg2_values;
    for (int index = 0; index <= static_cast<int>(std::acosh(std::sqrt(eps_r_prime)) / 0.01); ++index)
    {
        const double t = static_cast<double>(index) * 0.01;
        reg2_values.push_back(F_reg2(t, context));
    }
    double reg2_max = 1.0;
    for (const auto &value : reg2_values)
        reg2_max = std::max(reg2_max, std::abs(value));

    for (int index = 0; index <= static_cast<int>(std::acosh(std::sqrt(eps_r_prime)) / 0.01); ++index)
    {
        const double t = static_cast<double>(index) * 0.01;
        const double lambda = k_0 * std::cosh(t);
        const std::complex<double> f_full = F(lambda, context);
        const std::complex<double> f_sing =
            context.residue / (lambda - context.lambda_p);
        const std::complex<double> reg2_value = F_reg2(t, context);
        region2_t << t << "," << std::abs(f_full) << "," << std::abs(f_sing) << ","
                  << reg2_value.real() / reg2_max << "\n";
    }

    std::ofstream region3(out_dir / "scalar_pot_region3.csv");
    region3 << "lambda_over_k0,re_full,re_minus_static,im_full,im_minus_static\n";
    region3 << std::scientific << std::setprecision(10);
    const double region3_start_factor = std::sqrt(eps_r_prime) / 0.999;
    for (int index = 0; index <= static_cast<int>((10.0 - region3_start_factor) / 0.1); ++index)
    {
        const double lambda =
            (region3_start_factor + 0.1 * static_cast<double>(index)) * k_0;
        const std::complex<double> f_full = F(lambda, context);
        const std::complex<double> f_minus_static = F_reg3(lambda, context);
        region3 << lambda / k_0 << "," << f_full.real() << "," << f_minus_static.real()
                << "," << f_full.imag() << "," << f_minus_static.imag() << "\n";
    }

    std::ofstream summary(out_dir / "scalar_pot_summary.csv");
    summary << "key,value\n";
    summary << "freq_hz," << freq << "\n";
    summary << "lambda_0," << lambda_0 << "\n";
    summary << "eps_r_prime," << eps_r_prime << "\n";
    summary << "tan_delta," << tan_delta << "\n";
    summary << "k_0," << k_0 << "\n";
    summary << "k_real," << context.k.real() << "\n";
    summary << "k_imag," << context.k.imag() << "\n";
    summary << "d," << d << "\n";
    summary << "rho," << rho << "\n";
    summary << "lambda_p," << context.lambda_p << "\n";
    summary << "residue_real," << context.residue.real() << "\n";
    summary << "residue_imag," << context.residue.imag() << "\n";

    std::cout << "Arquivos gerados para scalar_pot em Cap_07/out/\n";
    return 0;
}
