// Cap_02/src/fdtd_1D_demo.cpp
#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <filesystem>
#include <fftw3.h>
#include <complex> // std::complex
#include <iomanip>





namespace fs = std::filesystem;

double compute_relative_error(const std::vector<std::complex<double>> &current,
                              const std::vector<std::complex<double>> &previous)
{
    double num = 0.0;   // Numerador: ||current - previous||^2
    double denom = 0.0; // Denominador: ||current||^2

    for (size_t i = 0; i < current.size(); ++i)
    {
        std::complex<double> diff = current[i] - previous[i];
        num += std::norm(diff);         // |diff|^2
        denom += std::norm(current[i]); // |current[i]|^2
    }

    if (denom == 0.0)
        return 0.0; // Evita divisão por zero (ou pode lançar exceção)

    return std::sqrt(num) / std::sqrt(denom);
}

int main()
{
    const std::string out_dir = PROJECT_OUT_DIR;
    std::cout << std::setprecision(15) << std::fixed;
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
    std::vector<std::vector<double>> V_period(M, std::vector<double>(Nz, 0.0));

    std::vector<std::vector<std::complex<double>>> V_prev_period_freq(M, std::vector<std::complex<double>>(Nz, {0.0, 0.0}));
    std::vector<std::vector<std::complex<double>>> V_period_freq(M, std::vector<std::complex<double>>(Nz, {0.0, 0.0}));
    std::cout << "1" << std::endl;
    // Time loop
    for (int nn = 2; nn <= Nk; ++nn)
    {
        std::cout << nn << std::endl;
        double Vo_nmin1 = V0 * std::cos(2.0 * M_PI * freq * (nn - 2.0) * delta_t); // Source.
        V_n[0] = (1.0 - beta1) * V_nmin1[0] - 2.0 * I_nmin1[0] + (2.0 / Rs) * Vo_nmin1;

        // Loop code, the vectorial instructions are not used in this traduction.
        for (int kk = 1; kk < (Nz - 1); ++kk)
            V_n[kk] = V_nmin1[kk] - (I_nmin1[kk] - I_nmin1[kk - 1]);

        V_n[Nz - 1] = (1.0 - beta2) * V_nmin1[Nz - 1] + 2.0 * I_nmin1[Nz - 2];

        // Loop code, the vectorial instructions are not used in this traduction.
        for (int kk = 0; kk < (Nz - 1); ++kk)
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
        int index = nn % M;
        if (index == 0)
            index = M;
        //std::cout << "index: " << index << std::endl;
        V_period[index] = V_n;

        //std::cout << "1.1" << std::endl;

        //std::cout << "V_period.size(): " << V_period.size() << std::endl;
        //std::cout << "V_period[0].size(): " << V_period[0].size() << std::endl;

        if (index == M)
        {
            // FFT over M time samples per spatial point (per column)
            // std::vector<std::complex<double>> V_period_freq(Nz, std::complex<double>(0.0, 0.0));
            // Aloca in/out e cria plano FFTW uma vez só
            fftw_complex *in = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * M);
            fftw_complex *out = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * M);
            fftw_plan p = fftw_plan_dft_1d(M, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
            //std::cout << "1.2" << std::endl;
            for (int z = 0; z < Nz; ++z)
            {
                std::cout << " = " << z << std::endl;
                // Preenche vetor de entrada com a coluna z
                for (int m = 0; m < M; ++m)
                {
                    in[m][0] = V_period[z][m]; // real
                    in[m][1] = 0.0;            // imag
                }
                //std::cout << "1.2.1" << std::endl;
                // Executa FFT
                fftw_execute(p);
                //std::cout << "1.2.2" << std::endl;
                for (int m = 0; m < M; ++m)
                {
                    V_period_freq[z][m] = std::complex<double>(out[m][0], out[m][1]); // real, imag
                }
                //std::cout << "1.2.3" << std::endl;
                // out[2][0]; // real part at k=2
                // out[2][1]; // imag part at k=2
            }
            //std::cout << "1.3" << std::endl;
            // Finaliza FFTW
            fftw_destroy_plan(p);
            fftw_free(in);
            fftw_free(out);
            //std::cout << "1.4" << std::endl;
            const double k = 1; // See discussion toward end of file regarding the index k.
            double eps = compute_relative_error(V_period_freq[k], V_prev_period_freq[k]);
            // Note that RMS norm includes inverse of root of length of vector, but it cancels above.
            // The FFTs in the numerator and denominator of the above expression are
            // both unscaled - the scale factors cancel here. See later comments regarding correct scaling of the FFT.

            // Exit loop, or overwrite for next period:
            if (eps < 0.002)
                break;
            //std::cout << "1.5" << std::endl;
            V_prev_period_freq = V_period_freq;
        }
    }
    //std::cout << "2" << std::endl;
    /* Now compute exact reults (p.36) and compare to the simulated ones.
    For the phasor results, what is needed is the first harmonic of the
    Fourier series expansion. The discussion in the textbook on p.58 refers to the general use of the
    FFT to approximate the Fourier transform. In this case, we are actually
    doing a Fourier series analysis. By sampling over exactly one period of the
    (eventually) periodic signal, the DFT (and hence FFT) exactly represent
    the (sampled) Fourier series. See, for instance, the discussion in E.O.Brigham, "The
    Fast Fourier Transform and its applications", Prentice-Hall 1988,
    p.98-100 where this is discussed in detail.
    In this case, there is an additional factor of the period 1/T (present in the Fourier
    series coefficients, but absent in the Fourier integral) which must be
    taken in to account.
    Furthermore, the factor of 2 comes from the negative and positive frequency
    components of the Fourier integral.*/
    //std::cout << "2" << std::endl;
    for (int i = 0; i < Nz; ++i)
    {
        V_n[i] *= delta_t / (C * delta_z);
        std::cout << " V_n[" << i << "]: " <<  V_n[i] << std::endl;
        for (int m = 0; m < M; ++m)
        {
            V_period_freq[m][i] *= delta_t / (C * delta_z);
            V_period_freq[m][i] *= (2.0 * delta_t / T);
        }
    }
    //std::cout << "3" << std::endl;
    const double k = 1;
    std::vector<double> z(Nz, 0.0);
    for (int i = 0; i < Nz; ++i)
    {
        z[i] = delta_z * i;
    }
    const double lambda = c / freq;
    const double beta = 2.0 * M_PI / lambda;
    const double Gamma = (Rl - Z_0) / (Rl + Z_0); // Eq. (2.16)
    const double V_plus = 0.5 * V0;               // for matched source, Eq. (2.15)

    std::vector<double> z_exact;
    for (double zz = 0.0; zz <= h + 1e-9; zz += 0.01)
    { // +1e-9 para evitar erro de ponto flutuante
        z_exact.push_back(zz);
    }

    std::vector<std::complex<double>> V_exact(z_exact.size());

    for (size_t i = 0; i < z_exact.size(); ++i)
    {
        double phase_shift = beta * (z_exact[i] - h);
        std::complex<double> term1 = std::exp(std::complex<double>(0, -phase_shift));
        std::complex<double> term2 = Gamma * std::exp(std::complex<double>(0, phase_shift));
        V_exact[i] = V_plus * (term1 + term2); // Eq. (2.14)
        std::cout << "V_exact[" << i << "]: " << V_exact[i] << std::endl;
    }

    std::ofstream file(out_dir + "/comparison_voltage.csv");

    // Exporta cabeçalho
    file << "z_exact,Re(V_exact),Im(V_exact),z_fdtd,Re(V_fdtd),Im(V_fdtd)\n";

    // Exporta z_exact e V_exact (pode ter tamanho diferente de Nz)
    size_t N_exact = z_exact.size();
    size_t N_fdtd = z.size(); // igual Nz

    size_t N_max = std::max(N_exact, N_fdtd);

    for (size_t i = 0; i < N_max; ++i)
    {
        if (i < N_exact)
            file << z_exact[i] << "," << real(V_exact[i]) << "," << imag(V_exact[i]);
        else
            file << ",,"; // Deixa vazio se não tem mais z_exact

        file << ",";

        if (i < N_fdtd)
            file << z[i] << "," << V_period_freq[k][i].real() << "," << V_period_freq[k][i].imag();
        else
            file << ",,";

        file << "\n";
    }

    file.close();

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
        file3 << V_period_freq[k][z].real() << "," << V_period_freq[k][z].imag() << "\n";
    file3.close();

    std::cout << "Simulation complete. Output saved in Cap_02/out/.\n";
    return 0;
}
