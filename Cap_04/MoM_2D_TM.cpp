#include "include/mom_tm_solver.hpp"
#include "include/cyl_tm_echo_width.hpp"

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numbers>
#include <stdexcept>
#include <string>
#include <vector>

namespace fs = std::filesystem;

namespace
{
// Translation of MoM_2D_TM.m:
// compute TM scattering from a PEC cylinder using the MoM EFIE with
// either a single-point off-diagonal evaluation or a simple quadrature.
struct Config
{
    std::string analysis_type = "current";
    int segments = 10;
    int n_per_lambda = 10;
    double lambda = 1.0;
    double radius = 1.0 / (2.0 * std::numbers::pi);
    double phi_inc = std::numbers::pi;
    double e0 = 1.0;
    int n_bessel = 25;
    double k_start_over_a = 0.1;
    double k_stop_over_a = 10.0;
    double k_step_over_a = 0.05;
    bool toeplitz = true;
};

std::string require_value(int &index, int argc, char **argv)
{
    if (index + 1 >= argc)
        throw std::runtime_error(std::string("Falta valor para ") + argv[index]);
    return argv[++index];
}

void print_help()
{
    std::cout
        << "Uso: ./build/MoM_2D_TM [opcoes]\n"
        << "  --analysis-type current|rcs|1|2   tipo de analise (padrao current)\n"
        << "  --segments N                  numero de segmentos para analise current (padrao 10)\n"
        << "  --n-per-lambda N              discretizacao por comprimento de onda no sweep RCS (padrao 10)\n"
        << "  --lambda VALOR                comprimento de onda no caso current (padrao 1.0)\n"
        << "  --radius VALOR                raio do cilindro (padrao lambda/(2*pi) para current, 0.03 para rcs)\n"
        << "  --phi-inc VALOR               angulo de incidencia em radianos (padrao pi)\n"
        << "  --e0 VALOR                    amplitude do campo incidente (padrao 1.0)\n"
        << "  --n-bessel N                  numero de termos para a solucao exata RCS (padrao 25)\n"
        << "  --k-start-over-a VALOR        inicio do sweep em ka/a (padrao 0.1)\n"
        << "  --k-stop-over-a VALOR         fim do sweep em ka/a (padrao 10.0)\n"
        << "  --k-step-over-a VALOR         passo do sweep em ka/a (padrao 0.05)\n"
        << "  --toeplitz 0|1                usa simetria Toeplitz (padrao 1)\n"
        << "  --help                        mostra esta ajuda\n";
}
} // namespace

int main(int argc, char **argv)
{
    const fs::path out_dir = PROJECT_OUT_DIR;
    fs::create_directories(out_dir);

    Config config;
    bool radius_explicit = false;
    try
    {
        for (int i = 1; i < argc; ++i)
        {
            const std::string arg = argv[i];
            if (arg == "--help" || arg == "-h")
            {
                print_help();
                return 0;
            }
            if (arg == "--analysis-type")
                config.analysis_type = require_value(i, argc, argv);
            else if (arg == "--segments")
                config.segments = std::stoi(require_value(i, argc, argv));
            else if (arg == "--n-per-lambda")
                config.n_per_lambda = std::stoi(require_value(i, argc, argv));
            else if (arg == "--lambda")
                config.lambda = std::stod(require_value(i, argc, argv));
            else if (arg == "--radius")
            {
                config.radius = std::stod(require_value(i, argc, argv));
                radius_explicit = true;
            }
            else if (arg == "--phi-inc")
                config.phi_inc = std::stod(require_value(i, argc, argv));
            else if (arg == "--e0")
                config.e0 = std::stod(require_value(i, argc, argv));
            else if (arg == "--n-bessel")
                config.n_bessel = std::stoi(require_value(i, argc, argv));
            else if (arg == "--k-start-over-a")
                config.k_start_over_a = std::stod(require_value(i, argc, argv));
            else if (arg == "--k-stop-over-a")
                config.k_stop_over_a = std::stod(require_value(i, argc, argv));
            else if (arg == "--k-step-over-a")
                config.k_step_over_a = std::stod(require_value(i, argc, argv));
            else if (arg == "--toeplitz")
                config.toeplitz = (std::stoi(require_value(i, argc, argv)) != 0);
            else
                throw std::runtime_error("Opcao desconhecida: " + arg);
        }
    }
    catch (const std::exception &ex)
    {
        std::cerr << "Erro: " << ex.what() << "\n";
        print_help();
        return 1;
    }

    const double eps_0 = 8.854e-12;
    const double mu_0 = 4.0 * std::numbers::pi * 1e-7;
    const double eta_0 = std::sqrt(mu_0 / eps_0);
    const double h_0 = config.e0 / eta_0;

    if (config.analysis_type == "1")
        config.analysis_type = "current";
    else if (config.analysis_type == "2")
        config.analysis_type = "rcs";

    if (config.analysis_type == "current")
    {
        if (!radius_explicit)
            config.radius = config.lambda / (2.0 * std::numbers::pi);

        const double k = 2.0 * std::numbers::pi / config.lambda;
        const MomTmResult current_mid =
            mom_tm_solver(k, config.segments, config.radius, config.e0, config.phi_inc, false, eta_0, config.toeplitz);
        const MomTmResult current_quad =
            mom_tm_solver(k, config.segments, config.radius, config.e0, config.phi_inc, true, eta_0, config.toeplitz);

        std::ofstream out(out_dir / "mom_tm_current.csv");
        out << "phi_deg,abs_js_single_over_h_inc,abs_js_quad_over_h_inc,"
               "re_js_single,im_js_single,re_js_quad,im_js_quad\n";
        out << std::scientific << std::setprecision(10);
        for (int i = 0; i < config.segments; ++i)
        {
            const double phi_deg = current_mid.phi_c[i] * 180.0 / std::numbers::pi;
            out << phi_deg << ","
                << std::abs(current_mid.current[i]) / h_0 << ","
                << std::abs(current_quad.current[i]) / h_0 << ","
                << current_mid.current[i].real() << ","
                << current_mid.current[i].imag() << ","
                << current_quad.current[i].real() << ","
                << current_quad.current[i].imag() << "\n";
        }

        std::ofstream summary(out_dir / "mom_tm_current_summary.csv");
        summary << "key,value\n";
        summary << "segments," << config.segments << "\n";
        summary << "lambda," << config.lambda << "\n";
        summary << "k," << k << "\n";
        summary << "radius," << config.radius << "\n";
        summary << "phi_inc," << config.phi_inc << "\n";
        summary << "e0," << config.e0 << "\n";
        summary << "toeplitz," << (config.toeplitz ? 1 : 0) << "\n";
        summary << "cond_single," << current_mid.condition_number << "\n";
        summary << "cond_quad," << current_quad.condition_number << "\n";

        std::cout << "Analise current concluida. Arquivos: mom_tm_current.csv e mom_tm_current_summary.csv\n";
    }
    else if (config.analysis_type == "rcs")
    {
        if (!radius_explicit)
            config.radius = 0.03;

        std::vector<double> k_values;
        for (double factor = config.k_start_over_a; factor <= config.k_stop_over_a + 1e-12; factor += config.k_step_over_a)
            k_values.push_back(factor / config.radius);

        const std::size_t n_freq = k_values.size();
        std::vector<int> segments_per_freq(n_freq, 0);
        std::vector<double> lambda(n_freq, 0.0);
        std::vector<double> rcs_tm(n_freq, 0.0);
        std::vector<double> rcs_tm_quad(n_freq, 0.0);
        std::vector<double> cond_single(n_freq, 0.0);
        std::vector<double> cond_quad(n_freq, 0.0);

        for (std::size_t ff = 0; ff < n_freq; ++ff)
        {
            std::cout << "Frequency " << (ff + 1) << " of " << n_freq << "\n";
            lambda[ff] = 2.0 * std::numbers::pi / k_values[ff];
            segments_per_freq[ff] = std::max(static_cast<int>(
                                                 std::ceil(2.0 * std::numbers::pi * config.radius / lambda[ff] *
                                                           config.n_per_lambda)),
                                             20);

            const MomTmResult single =
                mom_tm_solver(k_values[ff], segments_per_freq[ff], config.radius, config.e0,
                              config.phi_inc, false, eta_0, config.toeplitz);
            const MomTmResult quad =
                mom_tm_solver(k_values[ff], segments_per_freq[ff], config.radius, config.e0,
                              config.phi_inc, true, eta_0, config.toeplitz);
            cond_single[ff] = single.condition_number;
            cond_quad[ff] = quad.condition_number;

            std::complex<double> rcs_acc_single = 0.0;
            std::complex<double> rcs_acc_quad = 0.0;
            const double phi_bistat = config.phi_inc;
            for (int nn = 0; nn < segments_per_freq[ff]; ++nn)
            {
                const std::complex<double> phase =
                    std::exp(std::complex<double>(0.0, 1.0) *
                             k_values[ff] * (single.x_c[nn] * std::cos(phi_bistat) +
                                             single.y_c[nn] * std::sin(phi_bistat)));
                rcs_acc_single += single.current[nn] * single.strip_width * phase;
                rcs_acc_quad += quad.current[nn] * quad.strip_width * phase;
            }

            rcs_tm[ff] = (k_values[ff] * eta_0 * eta_0 / 4.0) * std::norm(rcs_acc_single) /
                         (config.e0 * config.e0);
            rcs_tm_quad[ff] = (k_values[ff] * eta_0 * eta_0 / 4.0) * std::norm(rcs_acc_quad) /
                              (config.e0 * config.e0);
        }

        const std::vector<double> exact =
            cyl_tm_echo_width(config.radius, k_values, config.n_bessel, std::numbers::pi);

        std::ofstream out(out_dir / "mom_tm_rcs.csv");
        out << "ka,a_over_lambda,rcs_quad_over_pi_a,rcs_single_over_pi_a,"
               "exact_over_lambda,rcs_quad_over_lambda,rcs_single_over_lambda,"
               "cond_single,cond_quad,segments\n";
        out << std::scientific << std::setprecision(10);
        for (std::size_t ff = 0; ff < n_freq; ++ff)
        {
            out << k_values[ff] * config.radius << ","
                << config.radius / lambda[ff] << ","
                << rcs_tm_quad[ff] / (std::numbers::pi * config.radius) << ","
                << rcs_tm[ff] / (std::numbers::pi * config.radius) << ","
                << exact[ff] / lambda[ff] << ","
                << rcs_tm_quad[ff] / lambda[ff] << ","
                << rcs_tm[ff] / lambda[ff] << ","
                << cond_single[ff] << ","
                << cond_quad[ff] << ","
                << segments_per_freq[ff] << "\n";
        }

        std::ofstream summary(out_dir / "mom_tm_rcs_summary.csv");
        summary << "key,value\n";
        summary << "radius," << config.radius << "\n";
        summary << "phi_inc," << config.phi_inc << "\n";
        summary << "e0," << config.e0 << "\n";
        summary << "n_per_lambda," << config.n_per_lambda << "\n";
        summary << "n_bessel," << config.n_bessel << "\n";
        summary << "k_start_over_a," << config.k_start_over_a << "\n";
        summary << "k_stop_over_a," << config.k_stop_over_a << "\n";
        summary << "k_step_over_a," << config.k_step_over_a << "\n";
        summary << "toeplitz," << (config.toeplitz ? 1 : 0) << "\n";

        std::cout << "Analise rcs concluida. Arquivos: mom_tm_rcs.csv e mom_tm_rcs_summary.csv\n";
    }
    else
    {
        std::cerr << "Erro: analysis_type deve ser 1, 2, current ou rcs.\n";
        return 1;
    }

    return 0;
}
