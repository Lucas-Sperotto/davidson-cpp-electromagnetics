// Cap_02/src/fdtd_1D_demo.cpp
#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <filesystem>
#include <fftw3.h>

namespace fs = std::filesystem;

int main() {
    const std::string out_dir = PROJECT_OUT_DIR;
    fs::create_directories(out_dir);

    const double h = 0.25;
    const double C = 1.0;
    const double L = 1.0;
    const double c = 1.0 / std::sqrt(L * C);
    const double Z_0 = std::sqrt(L / C);
    const double Rs = 1.0;

    double Rl;
    std::cout << "Load resistance? (Z_0 = 1 Ohm) Default: 2 -> ";
    std::string input;
    std::getline(std::cin, input);
    Rl = input.empty() ? 2.0 : std::stod(input);

    const double V0 = 1.0;
    const double freq = 4.0;
    const int Nz = 11;
    const double delta_z = h / (Nz - 1);

    const double T = 1.0 / freq;
    const int M = 64;
    const double delta_t = T / M;
    const int Nk = 16 * M;
    const double growth = 1.5;

    const double beta1 = 2.0 * delta_t / (Rs * C * delta_z);
    const double beta2 = 2.0 * delta_t / (Rl * C * delta_z);
    const double r = (delta_t * delta_t) / (L * C * delta_z * delta_z);

    std::vector<double> V_nmin1(Nz, 0.0), I_nmin1(Nz, 0.0);
    std::vector<double> V_n(Nz, 0.0), I_n(Nz, 0.0);
    std::vector<std::vector<double>> V_time_series;
    std::vector<std::vector<double>> V_period;

    for (int nn = 2; nn <= Nk; ++nn) {
        double Vo_nmin1 = V0 * std::cos(2 * M_PI * freq * (nn - 2) * delta_t);
        V_n[0] = (1 - beta1) * V_nmin1[0] - 2 * I_nmin1[0] + (2.0 / Rs) * Vo_nmin1;

        for (int kk = 1; kk < Nz - 1; ++kk)
            V_n[kk] = V_nmin1[kk] - (I_nmin1[kk] - I_nmin1[kk - 1]);

        V_n[Nz - 1] = (1 - beta2) * V_nmin1[Nz - 1] + 2 * I_nmin1[Nz - 2];

        for (int kk = 0; kk < Nz - 1; ++kk)
            I_n[kk] = I_nmin1[kk] - r * (V_n[kk + 1] - V_n[kk]);

        double norm_new = 0.0, norm_old = 0.0;
        for (int i = 0; i < Nz; ++i) {
            norm_new += V_n[i] * V_n[i];
            norm_old += V_nmin1[i] * V_nmin1[i];
        }

        if (nn > 5 && std::sqrt(norm_new) > growth * std::sqrt(norm_old)) {
            std::cerr << "Algorithm unstable. Time step: " << nn << "\n";
            break;
        }

        V_nmin1 = V_n;
        I_nmin1 = I_n;
        V_time_series.push_back(V_n);

        if (V_period.size() < M) {
            V_period.push_back(V_n);
        } else {
            V_period.erase(V_period.begin());
            V_period.push_back(V_n);
        }
    }

    // Export voltage at last time step
    std::ofstream file1(out_dir + "/fdtd_voltage.csv");
    for (const auto& v : V_n) file1 << v << "\n";
    file1.close();

    // Export voltage over time
    std::ofstream file2(out_dir + "/fdtd_time_series.csv");
    for (const auto& row : V_time_series) {
        for (size_t j = 0; j < row.size(); ++j)
            file2 << row[j] << (j == row.size() - 1 ? "\n" : ",");
    }
    file2.close();

    // FFT over M time samples per spatial point (per column)
    std::vector<std::vector<double>> spectrum(Nz, std::vector<double>(2, 0.0));
    for (int z = 0; z < Nz; ++z) {
        fftw_complex *in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * M);
        fftw_complex *out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * M);
        fftw_plan p = fftw_plan_dft_1d(M, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

        for (int m = 0; m < M; ++m) {
            in[m][0] = V_period[m][z]; // real
            in[m][1] = 0.0;            // imag
        }

        fftw_execute(p);
        spectrum[z][0] = out[2][0]; // real part at k=2
        spectrum[z][1] = out[2][1]; // imag part at k=2

        fftw_destroy_plan(p);
        fftw_free(in);
        fftw_free(out);
    }

    std::ofstream file3(out_dir + "/fdtd_spectrum.csv");
    for (int z = 0; z < Nz; ++z)
        file3 << spectrum[z][0] << "," << spectrum[z][1] << "\n";
    file3.close();

    std::cout << "Simulation complete. Output saved in Cap_02/out/.\n";
    return 0;
}
