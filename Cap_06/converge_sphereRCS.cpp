#include "sphereRCS.hpp"

#include "mom3d_globals.hpp"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

double require_double(const std::string& flag, const char* value)
{
    if (value == nullptr) {
        throw std::runtime_error("Missing value for " + flag);
    }
    return std::stod(value);
}

void print_help()
{
    std::cout
        << "Uso: ./build/converge_sphereRCS [opcoes]\n"
        << "  --a <m>          Raio da esfera (default: 1)\n"
        << "  --k-start <val>  Inicio de k (default: 0.1)\n"
        << "  --k-stop <val>   Fim de k (default: 10.0)\n"
        << "  --k-step <val>   Passo de k (default: 0.05)\n"
        << "  --help           Mostra esta ajuda\n";
}

}  // namespace

int main(int argc, char** argv)
{
    double a = 1.0;
    double k_start = 0.1;
    double k_stop = 10.0;
    double k_step = 0.05;

    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--a") {
            a = require_double(arg, (i + 1 < argc) ? argv[++i] : nullptr);
        } else if (arg == "--k-start") {
            k_start = require_double(arg, (i + 1 < argc) ? argv[++i] : nullptr);
        } else if (arg == "--k-stop") {
            k_stop = require_double(arg, (i + 1 < argc) ? argv[++i] : nullptr);
        } else if (arg == "--k-step") {
            k_step = require_double(arg, (i + 1 < argc) ? argv[++i] : nullptr);
        } else if (arg == "--help") {
            print_help();
            return 0;
        } else {
            throw std::runtime_error("Unknown argument: " + arg);
        }
    }

    if (k_step <= 0.0 || k_stop < k_start) {
        throw std::runtime_error("Invalid k range.");
    }

    std::filesystem::create_directories(cap06::project_out_dir());

    std::vector<double> k;
    for (double value = k_start; value <= k_stop + 1.0e-12; value += k_step) {
        k.push_back(value);
    }

    const std::vector<double> rcs_1 = cap06::sphereRCS(a, k, 1);
    const std::vector<double> rcs_5 = cap06::sphereRCS(a, k, 5);
    const std::vector<double> rcs_10 = cap06::sphereRCS(a, k, 10);
    const std::vector<double> rcs_50 = cap06::sphereRCS(a, k, 50);

    std::ofstream curve_file(cap06::project_out_dir() / "sphere_rcs_convergence.csv");
    curve_file << "k,ka,a_over_lambda,rcs_n1_over_pi_a2,rcs_n5_over_pi_a2,rcs_n10_over_pi_a2,rcs_n50_over_pi_a2\n";
    for (std::size_t i = 0; i < k.size(); ++i) {
        const double ka = k[i] * a;
        const double lambda = 2.0 * cap06::PI / k[i];
        curve_file << k[i] << ',' << ka << ',' << a / lambda << ','
                   << rcs_1[i] / (cap06::PI * a * a) << ','
                   << rcs_5[i] / (cap06::PI * a * a) << ','
                   << rcs_10[i] / (cap06::PI * a * a) << ','
                   << rcs_50[i] / (cap06::PI * a * a) << '\n';
    }

    std::ofstream summary_file(cap06::project_out_dir() / "sphere_rcs_summary.csv");
    summary_file << "key,value\n";
    summary_file << "a_m," << a << '\n';
    summary_file << "k_start," << k_start << '\n';
    summary_file << "k_stop," << k_stop << '\n';
    summary_file << "k_step," << k_step << '\n';
    summary_file << "num_samples," << k.size() << '\n';

    std::cout << "converge_sphereRCS finalizado com " << k.size() << " amostras.\n";
    return 0;
}
