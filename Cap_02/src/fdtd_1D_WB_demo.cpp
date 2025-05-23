// Cap_02/src/fdtd_1D_WB_demo.cpp

#include <iostream>   // Para entrada e saída padrão (ex: std::cout)
                      // Doc: https://en.cppreference.com/w/cpp/io
#include <vector>     // Para usar o container std::vector (vetores dinâmicos)
                      // Doc: https://en.cppreference.com/w/cpp/container/vector
#include <cmath>      // Funções matemáticas padrão (ex: std::sqrt, std::norm)
                      // Doc: https://en.cppreference.com/w/cpp/header/cmath
#include <fstream>    // Para operações de arquivo (ex: std::ofstream)
                      // Doc: https://en.cppreference.com/w/cpp/header/fstream
#include <filesystem> // Para criação e manipulação de diretórios
                      // Doc: https://en.cppreference.com/w/cpp/filesystem
#include <fftw3.h>    // Biblioteca FFTW para transformada rápida de Fourier
                      // Doc: http://www.fftw.org/doc/
#include <complex>    // Para números complexos com std::complex
                      // Doc: https://en.cppreference.com/w/cpp/numeric/complex
#include <iomanip>    // Para manipular formatação de saída (ex: std::setprecision)
                      // Doc: https://en.cppreference.com/w/cpp/io/manip

#include <string>

// Definindo um alias para a biblioteca filesystem para simplificar chamadas
namespace fs = std::filesystem;

// Função impulso gaussiano derivado normalizado
double gaussder_norm(double, double, double);

// Função principal do programa
int main()
{

    // Define o diretório de saída onde os arquivos de resultado serão salvos
    // O diretório de saída é definido como "Cap_02/out/"
    // essa informação é deninida no CMakeLists.txt
    // e é passada como variável de ambiente para o programa
    const std::string out_dir = PROJECT_OUT_DIR;

    // Define a precisão de saída padrão para 15 casas decimais e formato fixo (não científico)
    std::cout << std::setprecision(15) << std::fixed;

    // Cria o diretório de saída (se não existir)
    fs::create_directories(out_dir);

    std::ofstream voltage_log(out_dir + "/WB_voltage_over_time.csv");

    const double C = 1.0;                    // Capacitância por unidade de comprimento [F/m]
    const double L = 1.0;                    // Indutância por unidade de comprimento [H/m]
    const double Z_0 = std::sqrt(L / C);     // Impedância característica da linha de transmissão eq. (2.87)
    const double c = 1.0 / std::sqrt(L * C); // Velocidade de propagação na linha (v = 1/sqrt(LC))

    std::string resp;
    std::cout << "Análise por pulso gaussiano derivado (1) ou função transferência (0)? [1]: ";
    std::getline(std::cin, resp);
    bool flag = resp.empty() || resp == "1";

    double h, Rs, Rl, Courant, sigma;
    int Nk;

    // For time-domain demo, use Rs=1, Rl=3 for nice results. Start out with both %=1
    if (flag)
    {
        h = 0.5;  // Comprimento total da linha de transmissão (em metros) [m]
        Rs = 1.0; // Resistência da fonte [Ohm]
        std::cout << "Resistência de carga (default 3 Ohm): ";
        std::getline(std::cin, resp);
        Rl = resp.empty() ? 3.0 : std::stod(resp);
        Courant = 0.5; // 0.7
        sigma = 1.0 / (8.0 * M_PI);
        Nk = 512; // Number of time steps. Power of two best (for FFT).
    }
    else
    {
        h = 0.25; // Comprimento total da linha de transmissão (em metros) [m]
        Rs = 0.5; // Resistência da fonte [Ohm]
        Rl = 2.0;
        Courant = 0.4;
        sigma = 1.0 / (16.0 * M_PI); // more bandwidth
        Nk = 1024;                   // Number of time steps. Power of two best (for FFT).
    }

    const int Nz = 11;                   // 41; % 11, 21 and 41 demonstrate effect of dispersion nicely.
    const double delta_z = h / (Nz - 1); // Number of points
    const double delta_t = Courant * delta_z / c;
    const double offset = 4.0 * sigma;

    const double beta1 = 2.0 * delta_t / (Rs * C * delta_z); // Eq. (2.67)
    const double beta2 = 2.0 * delta_t / (Rl * C * delta_z); // Eq. (2.68)
    const double const_beta = delta_t / (C * delta_z);
    // const double r = (delta_t * delta_t) / (L * C * delta_z * delta_z); // Eq. (2.69)

    const double growth = 20.0; // An abritary growth factor indicating instability
    // Don't make too close to 1, since the early time
    // behaviour can indeed grow quite quickly.

    std::vector<double> V_nmin1(Nz, 0.0), I_nmin1(Nz, 0.0);
    std::vector<double> V_n(Nz, 0.0), I_n(Nz, 0.0); // Eq. (2.61) e Eq. (2.62)
    std::vector<double> Vs(Nk, 0.0), Vl(Nk, 0.0), time(Nk, 0.0);

    // Cabeçalho: TimeStep, Vn[1], Vn[2], ..., Vn[Nz]
    voltage_log << "TimeStep";
    for (size_t i = 0; i < Nz; ++i)
        voltage_log << ",Vn[" << (i + 1) << "]";
    voltage_log << "\n";

    for (int nn = 1; nn < Nk; ++nn)
    {
        double t = (nn)*delta_t;
        time[nn] = t;
    }

    for (int nn = 1; nn < Nk; ++nn)
    {
        double Vo_nmin1 = gaussder_norm(time[nn], offset, sigma);
        Vs[nn] = Vo_nmin1;

        V_n[0] = (1.0 - beta1) * V_nmin1[0] - 2.0 * const_beta * I_nmin1[0] + (2.0 * const_beta / Rs) * Vo_nmin1;

        for (int kk = 1; kk < Nz - 1; ++kk)
            V_n[kk] = V_nmin1[kk] - const_beta * (I_nmin1[kk] - I_nmin1[kk - 1]);

        V_n[Nz - 1] = (1.0 - beta2) * V_nmin1[Nz - 1] + 2.0 * const_beta * I_nmin1[Nz - 2];
        Vl[nn] = V_n[Nz - 1];

        for (int kk = 0; kk < Nz - 1; ++kk)
            I_n[kk] = I_nmin1[kk] - delta_t / (L * delta_z) * (V_n[kk + 1] - V_n[kk]);

        voltage_log << nn;
        for (size_t i = 0; i < Nz; ++i)
        {
            double voltage = (delta_t / (C * delta_z)) * V_n[i];
            voltage_log << "," << voltage;
        }
        voltage_log << "\n";

        double norm_new = 0.0, norm_old = 0.0;
        for (int i = 0; i < Nz; ++i)
        {
            norm_new += V_n[i] * V_n[i];
            norm_old += V_nmin1[i] * V_nmin1[i];
        }
        if (nn > 10 && std::sqrt(norm_new) > growth * std::sqrt(norm_old))
        {
            std::cerr << "Instabilidade detectada em nn = " << nn << "\n";
            break;
        }
        V_nmin1 = V_n;
        I_nmin1 = I_n;
    }

    double tau_s = Z_0 / (Rs + Z_0);

    voltage_log.close();

    // for (int i = 0; i < Nk; ++i)
    // std::cout << "time[" << i << "]" << time[i] << std::endl;
    // for (int i = 0; i < Nk; ++i)
    // std::cout << "Vs[" << i << "]" << Vs[i] << std::endl;
    // for (int i = 0; i < Nz; ++i)
    // std::cout << "V_n[" << i << "]" << V_n[i] << std::endl;
    // for (int i = 0; i < Nk; ++i)
    // std::cout << "Vl[" << i << "]" << Vl[i] << std::endl;

    if (flag)
    {
        std::ofstream f_vs(out_dir + "/fdtd_wb_vs.csv"), f_vl(out_dir + "/fdtd_wb_vl.csv");

        for (int i = 0; i < Nk; ++i)
        {
            f_vs << time[i] << "," << Vs[i] * tau_s << "\n";
            f_vl << time[i] << "," << Vl[i] << "\n";
        }
        f_vs.close();
        f_vl.close();
    }
    else
    {
        fftw_complex *in1 = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * Nk);
        fftw_complex *in2 = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * Nk);
        fftw_complex *out1 = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * Nk);
        fftw_complex *out2 = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * Nk);

        fftw_plan p1 = fftw_plan_dft_1d(Nk, in1, out1, FFTW_FORWARD, FFTW_ESTIMATE);
        fftw_plan p2 = fftw_plan_dft_1d(Nk, in2, out2, FFTW_FORWARD, FFTW_ESTIMATE);

        for (int i = 0; i < Nk; ++i)
        {
            in1[i][0] = Vl[i];
            in1[i][1] = 0;
            in2[i][0] = Vs[i];
            in2[i][1] = 0;
        }
        fftw_execute(p1);
        fftw_execute(p2);

        std::ofstream tf_file(out_dir + "/fdtd_wb_tf.csv"), tf_file_vl(out_dir + "/fdtd_wb_tf_vl.csv");

        for (int i = 0; i < Nk; ++i)
        {
            tf_file_vl << time[i] << "," << Vl[i] << "\n";
        }
        tf_file_vl.close();
        std::vector<std::complex<double>> Tx(Nk, 0.0);
        double delta_f = 1.0 / (Nk * delta_t);
        for (int i = 0; i < (Nk / 32); ++i)
        {
            double f = i * delta_f;
            std::complex<double> G(out1[i][0], out1[i][1]);
            std::complex<double> H(out2[i][0], out2[i][1]);
            Tx[i] = G / H;
            // std::cout << "Tx[" << i << "]" << Tx[i] << std::endl;
            tf_file << f << "," << std::abs(Tx[i]) << "\n";
        }
        tf_file.close();
        fftw_destroy_plan(p1);
        fftw_destroy_plan(p2);
        fftw_free(in1);
        fftw_free(in2);
        fftw_free(out1);
        fftw_free(out2);
    }
    /*
        // std::cout << "ans: " << ans << std::endl;
        std::cout << "beta1: " << beta1 << std::endl;
        std::cout << "beta2: " << beta2 << std::endl;
        std::cout << "c: " << c << std::endl;
        std::cout << "C: " << C << std::endl;
        std::cout << "const_beta: " << const_beta << std::endl;
        std::cout << "Courant: " << Courant << std::endl;
        std::cout << "delta_t: " << delta_t << std::endl;
        std::cout << "delta_z: " << delta_z << std::endl;
        std::cout << "flag: " << flag << std::endl;
        std::cout << "growth: " << growth << std::endl;
        std::cout << "h: " << h << std::endl;
        std::cout << "L: " << L << std::endl;
        std::cout << "ofset: " << offset << std::endl;
        std::cout << "sigma: " << sigma << std::endl;
        // std::cout << "r: " << r << std::endl;
        std::cout << "Rl: " << Rl << std::endl;
        std::cout << "Rs: " << Rs << std::endl;
        std::cout << "tau_s: " << tau_s << std::endl;
        std::cout << "Z_0: " << Z_0 << std::endl;
    */
    std::cout << "Simulação banda larga finalizada. Resultados salvos em Cap_02/out/.\n";
    return 0;
}

// Função impulso gaussiano derivado normalizado
double gaussder_norm(double t, double m, double sigma)
{
    // return -1.0 / std::sqrt(2.0 * M_PI) * (t - m) / sigma * sigma * sigma * exp(-(t-m)^2/(2*sigma^2));
    return -std::exp(0.5) * (t - m) / sigma * std::exp(-std::pow((t - m), 2) / (2 * sigma * sigma));
}