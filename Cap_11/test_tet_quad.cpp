#include "include/tet_quad.hpp"

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

namespace fs = std::filesystem;

namespace
{
std::string require_value(int &index, int argc, char **argv)
{
    if (index + 1 >= argc)
        throw std::runtime_error(std::string("Falta valor para ") + argv[index]);
    return argv[++index];
}
}

int main(int argc, char **argv)
{
    int n = 11;
    for (int i = 1; i < argc; ++i)
    {
        const std::string arg = argv[i];
        if (arg == "--n")
            n = std::stoi(require_value(i, argc, argv));
    }

    // Translation of test_tet_quad.m:
    // test the tetrahedral quadrature rule on low-order polynomials.
    const auto [w, lambda] = tet_quad(n);

    double quad_int = 0.0;
    for (int ii = 0; ii < n; ++ii)
        quad_int += w(ii) * lambda(ii, 0) * lambda(ii, 0) *
                    lambda(ii, 1) * lambda(ii, 1);
    const double quad_1 = quad_int;
    const double int_1 = 1.0 / 210.0;
    const double err_1 = std::abs(int_1 - quad_1) / quad_1;

    quad_int = 0.0;
    for (int ii = 0; ii < n; ++ii)
        quad_int += w(ii) * lambda(ii, 0) * lambda(ii, 1) *
                    lambda(ii, 2) * lambda(ii, 3);
    const double quad_2 = quad_int;
    const double int_2 = 1.0 / 840.0;
    const double err_2 = std::abs(int_2 - quad_2) / quad_2;

    quad_int = 0.0;
    for (int ii = 0; ii < n; ++ii)
        quad_int += w(ii);
    const double quad_3 = quad_int;
    const double int_3 = 1.0;
    const double err_3 = std::abs(int_3 - quad_3) / quad_3;

    const fs::path out_dir = PROJECT_OUT_DIR;
    fs::create_directories(out_dir);
    std::ofstream summary(out_dir / "test_tet_quad_summary.csv");
    summary << "test,reference,quadrature,relative_error\n";
    summary << std::scientific << std::setprecision(10);
    summary << "lambda1_sq_lambda2_sq," << int_1 << "," << quad_1 << "," << err_1 << "\n";
    summary << "lambda1_lambda2_lambda3_lambda4," << int_2 << "," << quad_2 << "," << err_2 << "\n";
    summary << "constant," << int_3 << "," << quad_3 << "," << err_3 << "\n";

    std::cout << "Arquivo gerado: test_tet_quad_summary.csv\n";
    return 0;
}
