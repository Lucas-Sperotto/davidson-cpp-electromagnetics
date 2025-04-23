// Cap_02/src/fdtd_1D_demo.cpp
#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <filesystem>
#include <fftw3.h>

namespace fs = std::filesystem;

// Calcula norma (equivalente MATLAB norm)
double fftw_norm(fftw_complex* x, int N) {
    double soma = 0.0;
    for (int k = 0; k < N; ++k) {
        double real = x[k][0];
        double imag = x[k][1];
        soma += real * real + imag * imag;
    }
    return sqrt(soma);
}




int main()
{
    const std::string out_dir = PROJECT_OUT_DIR;
    fs::create_directories(out_dir);

    const double h = 0.25; // length [m]
    const double C = 1.0;
    const double L = 1.0;
    const double c = 1.0 / std::sqrt(L * C);
    const double Z_0 = std::sqrt(L / C);
    const double Rs = 1.0; //[Ohm]

    double Rl;
    std::cout << "Load resistance? (Z_0 = 1 Ohm) Default: 2 -> ";
    std::string input;
    std::getline(std::cin, input);
    Rl = input.empty() ? 2.0 : std::stod(input);

    const double V0 = 1.0; // Amplitude of source voltage [V]
    // const double Rl = 1;
    const double freq = 4.0; // freq of applied sinusoid [Hz]
    const int Nz = 11;
    const double delta_z = h / (Nz - 1); // Number of points

    // Nk = 400; // Number of time steps for 0.25

    const double T = 1.0 / freq;  // Period of signal [s]
    const double delta_f = 1 / T; // See p.58 for discussion of FFT parameters.

    const int M = 64; // Samples per period.

    // const double delta_t = T / (M - 1);
    /* Note slight change. The reason is that M in this case is the
        number of points PER PERIOD, and we need to ensure that the last time sample in one period
        is not also the same time sample in the next period.*/
    const double delta_t = T / M;

    const int Nk = 16 * M; // A maximum number of periods to run if the convergence criteria eps is not achieved.

    /* An abritary growth factor indicating instability
        Don't make too close to 1, since the early time
        behaviour can indeed grow quite quickly.*/
    const double growth = 1.5;

    const double beta1 = 2.0 * delta_t / (Rs * C * delta_z);
    const double beta2 = 2.0 * delta_t / (Rl * C * delta_z);
    const double r = (delta_t * delta_t) / (L * C * delta_z * delta_z);

    //  First time step - Initialize.
    std::vector<double> V_nmin1(Nz, 0.0), I_nmin1(Nz, 0.0);
    // Pre-allocation
    std::vector<double> V_n(Nz, 0.0), I_n(Nz, 0.0);

    std::vector<std::vector<double>> V_time_series;
    std::vector<std::vector<double>> V_period;

    // Time loop
    for (int nn = 2; nn <= Nk; ++nn)
    {
        double Vo_nmin1 = V0 * std::cos(2.0 * M_PI * freq * (nn - 2.0) * delta_t); // Source.
        V_n[0] = (1.0 - beta1) * V_nmin1[0] - 2.0 * I_nmin1[0] + (2.0 / Rs) * Vo_nmin1;

        // Loop code, the vectorial instructions are not used in this traduction.
        for (int kk = 1; kk < Nz - 1; ++kk)
            V_n[kk] = V_nmin1[kk] - (I_nmin1[kk] - I_nmin1[kk - 1]);

        V_n[Nz - 1] = (1.0 - beta2) * V_nmin1[Nz - 1] + 2.0 * I_nmin1[Nz - 2];

        // Loop code, the vectorial instructions are not used in this traduction.
        for (int kk = 0; kk < Nz - 1; ++kk)
            I_n[kk] = I_nmin1[kk] - r * (V_n[kk + 1] - V_n[kk]);

        double norm_new = 0.0, norm_old = 0.0;
        for (int i = 0; i < Nz; ++i)
        {
            norm_new += V_n[i] * V_n[i];
            norm_old += V_nmin1[i] * V_nmin1[i];
        }

        if (nn > 5 && std::sqrt(norm_new) > growth * std::sqrt(norm_old))
        {
            std::cerr << "Algorithm unstable. Time step: " << nn << "\n";
            break;
        }

        V_nmin1 = V_n;
        I_nmin1 = I_n;
        V_time_series.push_back(V_n);

        if (V_period.size() < M)
        {
            V_period.push_back(V_n);
        }
        else
        {
            V_period.erase(V_period.begin());
            V_period.push_back(V_n);
        }

        // FFT over M time samples per spatial point (per column)
        std::vector<std::vector<double>> spectrum(Nz, std::vector<double>(2, 0.0));
        for (int z = 0; z < Nz; ++z)
        {
            fftw_complex *in = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * M);
            fftw_complex *out = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * M);
            fftw_plan p = fftw_plan_dft_1d(M, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

            for (int m = 0; m < M; ++m)
            {
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

        const double k = 2; / See discussion toward end of file regarding the index k.
        double eps = norm(V_period_freq(k,:) - V_prev_period_freq(k,:))/norm(V_period_freq(k,:))
        // Note that RMS norm includes inverse of root of length of vector, but it cancels above.  
        // The FFTs in the numerator and denominator of the above expression are
        // both unscaled - the scale factors cancel here. See later comments regarding correct scaling of the FFT.
        
        // Exit loop, or overwrite for next period:
        if (eps < 0.002) 
            break;

        V_prev_period_freq = V_period_freq;
    }

    // Export voltage at last time step
    std::ofstream file1(out_dir + "/fdtd_voltage.csv");
    for (const auto &v : V_n)
        file1 << v << "\n";
    file1.close();

    // Export voltage over time
    std::ofstream file2(out_dir + "/fdtd_time_series.csv");
    for (const auto &row : V_time_series)
    {
        for (size_t j = 0; j < row.size(); ++j)
            file2 << row[j] << (j == row.size() - 1 ? "\n" : ",");
    }
    file2.close();

    std::ofstream file3(out_dir + "/fdtd_spectrum.csv");
    for (int z = 0; z < Nz; ++z)
        file3 << spectrum[z][0] << "," << spectrum[z][1] << "\n";
    file3.close();

    std::cout << "Simulation complete. Output saved in Cap_02/out/.\n";
    return 0;
}
