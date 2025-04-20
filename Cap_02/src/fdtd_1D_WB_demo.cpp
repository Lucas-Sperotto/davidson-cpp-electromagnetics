// Cap_02/src/fdtd_1D_WB_demo.cpp
#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <string>
#include <complex>
#include <fftw3.h>
#include <filesystem>

namespace fs = std::filesystem;

// Função impulso gaussiano derivado normalizado
double gaussder_norm(double t, double m, double sigma) {
    return -std::exp(0.5) * (t - m) / sigma * std::exp(-std::pow((t - m), 2) / (2 * sigma * sigma));
}

int main() {
    const std::string out_dir = PROJECT_OUT_DIR;
    fs::create_directories(out_dir);

    const double C = 1.0;
    const double L = 1.0;
    const double Z_0 = std::sqrt(L / C);
    const double c = 1.0 / std::sqrt(L * C);

    std::string resp;
    std::cout << "Análise por pulso gaussiano derivado (1) ou função transferência (0)? [1]: ";
    std::getline(std::cin, resp);
    bool flag = resp.empty() || resp == "1";

    double h, Rs, Rl, Courant, sigma;
    int Nk;
    if (flag) {
        h = 0.5;
        Rs = 1.0;
        std::cout << "Resistência de carga (default 3 Ohm): ";
        std::getline(std::cin, resp);
        Rl = resp.empty() ? 3.0 : std::stod(resp);
        Courant = 0.5;
        sigma = 1.0 / (8.0 * M_PI);
        Nk = 512;
    } else {
        h = 0.25;
        Rs = 0.5;
        Rl = 2.0;
        Courant = 0.4;
        sigma = 1.0 / (16.0 * M_PI);
        Nk = 1024;
    }

    const int Nz = 11;
    const double delta_z = h / (Nz - 1);
    const double delta_t = Courant * delta_z / c;
    const double offset = 4.0 * sigma;

    const double beta1 = 2.0 * delta_t / (Rs * C * delta_z);
    const double beta2 = 2.0 * delta_t / (Rl * C * delta_z);
    const double const_beta = delta_t / (C * delta_z);
    const double growth = 20.0;

    std::vector<double> V_nmin1(Nz, 0.0), I_nmin1(Nz, 0.0);
    std::vector<double> V_n(Nz, 0.0), I_n(Nz, 0.0);
    std::vector<double> Vs(Nk, 0.0), Vl(Nk, 0.0), time(Nk, 0.0);

    for (int nn = 1; nn < Nk; ++nn) {
        double t = nn * delta_t;
        time[nn] = t;
        double Vo = gaussder_norm(t, offset, sigma);
        Vs[nn] = Vo;

        V_n[0] = (1 - beta1) * V_nmin1[0] - 2 * const_beta * I_nmin1[0] + (2 * const_beta / Rs) * Vo;

        for (int kk = 1; kk < Nz - 1; ++kk)
            V_n[kk] = V_nmin1[kk] - const_beta * (I_nmin1[kk] - I_nmin1[kk - 1]);

        V_n[Nz - 1] = (1 - beta2) * V_nmin1[Nz - 1] + 2 * const_beta * I_nmin1[Nz - 2];
        Vl[nn] = V_n[Nz - 1];

        for (int kk = 0; kk < Nz - 1; ++kk)
            I_n[kk] = I_nmin1[kk] - delta_t / (L * delta_z) * (V_n[kk + 1] - V_n[kk]);

        double norm_new = 0.0, norm_old = 0.0;
        for (int i = 0; i < Nz; ++i) {
            norm_new += V_n[i] * V_n[i];
            norm_old += V_nmin1[i] * V_nmin1[i];
        }
        if (nn > 10 && std::sqrt(norm_new) > growth * std::sqrt(norm_old)) {
            std::cerr << "Instabilidade detectada em nn = " << nn << "\n";
            break;
        }
        V_nmin1 = V_n;
        I_nmin1 = I_n;
    }

    std::ofstream f_vs(out_dir + "/fdtd_wb_vs.csv"), f_vl(out_dir + "/fdtd_wb_vl.csv");
    for (int i = 0; i < Nk; ++i) {
        f_vs << time[i] << "," << Vs[i] << "\n";
        f_vl << time[i] << "," << Vl[i] << "\n";
    }
    f_vs.close();
    f_vl.close();

    if (!flag) {
        fftw_complex *in1 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nk);
        fftw_complex *in2 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nk);
        fftw_complex *out1 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nk);
        fftw_complex *out2 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nk);

        fftw_plan p1 = fftw_plan_dft_1d(Nk, in1, out1, FFTW_FORWARD, FFTW_ESTIMATE);
        fftw_plan p2 = fftw_plan_dft_1d(Nk, in2, out2, FFTW_FORWARD, FFTW_ESTIMATE);

        for (int i = 0; i < Nk; ++i) {
            in1[i][0] = Vs[i]; in1[i][1] = 0;
            in2[i][0] = Vl[i]; in2[i][1] = 0;
        }
        fftw_execute(p1);
        fftw_execute(p2);

        std::ofstream tf_file(out_dir + "/fdtd_wb_tf.csv");
        double delta_f = 1.0 / (Nk * delta_t);
        for (int i = 0; i < Nk / 2; ++i) {
            double f = i * delta_f;
            std::complex<double> H(out2[i][0], out2[i][1]);
            std::complex<double> G(out1[i][0], out1[i][1]);
            std::complex<double> Tx = H / G;
            tf_file << f << "," << std::abs(Tx) << "\n";
        }
        tf_file.close();
        fftw_destroy_plan(p1);
        fftw_destroy_plan(p2);
        fftw_free(in1); fftw_free(in2);
        fftw_free(out1); fftw_free(out2);
    }

    std::cout << "Simulação banda larga finalizada. Resultados salvos em Cap_02/out/.\n";
    return 0;
}
