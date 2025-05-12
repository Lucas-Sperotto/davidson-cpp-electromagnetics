// src/fdtd2d.cpp
#include "../include/fdtd2d.hpp"
#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <filesystem>
#include "gaussder.cpp"

// Constantes físicas
const double c = 2.99792458E8;                //[m/s] Speed of light in vacuum
const double eps_0 = 8.854187817E-12;         //[F/m] epsilon_0
const double mu_0 = 4.0 * M_PI * 1E-7;        //[H/m] mu_0
const double eta_0 = std::sqrt(mu_0 / eps_0); //[Ohm] Wave impedance of free space

int run_fdtd2d_simulation()
{
    const std::string out_dir = PROJECT_OUT_DIR;

    // Define a precisão de saída padrão para 15 casas decimais e formato fixo (não científico)
    std::cout << std::setprecision(15) << std::fixed;

    // Cria o diretório de saída (se não existir)
    std::filesystem::create_directories(out_dir);


bool cyl_present = false; // Set to true if you want to include a PEC cylinder

    int refine = 1;
    std::cout << "Factor to refine mesh? 1 standard: ";
    std::cin >> refine;

    double pulse_compress = 1; // 2 % A quick way to shorten the pulse. For Figs 3.11-14, use 2. For 3.15, use 1

    int N_x = refine * 200;                             // refine*400     % number of cells in x-direction
    int N_y = refine * 100;                             // refine*200            % ditto y.
    int M = refine * 512;                               // refine*1024      % Number of time steps
    int L = std::round(static_cast<double>(N_x) / 2.0); // scat/tot field boundary on left side of Fig. 3.1 as Section 3.2
    int L1 = 20;                                        // 10              % New scat/tot field boundaries on upper, lower and right side of Fig. 3.1

    double delta_s = 0.005 * 2.0 / refine;               //[m] spatial step. This generates a physically larger computational volume.
    double R = 1.0;                                      // fraction of Courant limit. Must be <= 1
    double delta_t = R * delta_s / (c * std::sqrt(2.0)); //[s] Time step size
    double sigma = 1.0E-10 / pulse_compress;             // Controls spectral content of Gaussian derivative pulse equals 1/omega_max

    double f_max = (1.0 / sigma) / (2.0 * M_PI); // Freq of largest significant spectral component, in Hz. Information only
    double m_offset = 4.0 * sigma;               // Controls switch-on time
    double Peak = 1.0;                           // Peak amplitude of E field

    // Set parameters for PEC cylinder
    double radius = 0.03; // [m] radius of cylinder
    int N_centre_x = std::round(0.75 * N_x);
    int N_centre_y = std::round(0.5 * N_y);

    std::vector<double> E_y_point1(M);
    int point1_x = N_x / 4;
    int point1_y = N_y / 2;

    // Check that the simulation specification is valid:
    if ((N_centre_x - L) * delta_s <= radius)
    {
        std::cout << "Error in simulation data. Scattered/total field not entirely to the left of target" << std::endl;
        return 0;
    }

    // Set up material grid (free space to start)
    std::vector<std::vector<double>> C_Ex(N_x + 1, std::vector<double>(N_y + 1, delta_t / (eps_0 * delta_s)));
    std::vector<std::vector<double>> C_Ey(N_x + 1, std::vector<double>(N_y + 1, delta_t / (eps_0 * delta_s)));
    double D_Hz = delta_t / (mu_0 * delta_s);
    // Now force the electric fields to zero inside (and on the surface of) the PEC
    // Note that the indices of the centre are treated as per usual FDTD
    // indices, i.e. the actual location is:
    // x_c=(N_centre_x-1))*delta_s ; y_c=(N_centre_x-1))*delta_s
    if (cyl_present) // Otherwise just leave it as free space
    {
        for (int ii = 0; ii <= N_x; ++ii)
        {
            for (int jj = 0; jj <= N_y; ++jj)
            {
                if (sqrt(((ii - 0.5 - (N_centre_x - 1)) * delta_s) * ((ii - 0.5 - (N_centre_x - 1)) * delta_s) + ((jj - 1 - (N_centre_y - 1)) * delta_s) * ((jj - 1 - (N_centre_y - 1)) * delta_s)) <= radius)
                    C_Ex[ii][jj] = 0.0;

                if (sqrt(((ii - 1.0 - (N_centre_x - 1)) * delta_s) * ((ii - 1.0 - (N_centre_x - 1)) * delta_s) + ((jj - 0.5 - (N_centre_y - 1)) * delta_s) * ((jj - 0.5 - (N_centre_y - 1)) * delta_s)) <= radius)
                    C_Ey[ii][jj] = 0.0;
            }
        }
    }

    // Set up storage for time histories.
    std::vector<double> H_z_point1(M, 0.0);
    std::vector<double> E_y_point1(M, 0.0);
    int point1_x = N_x / 4;
    // This is another hack! also remove!!
    // point1_x = 100;
    int point1_y = N_y / 2;
    std::vector<double> H_z_point2(M, 0.0);
    std::vector<double> E_y_point2(M, 0.0);
    int point2_x = std::round(static_cast<double>(N_x + L) / 2.0);
    int point2_y = N_y / 2;

    // Produce a simple graphical output, showing the cylinder, scat/tot zone
    // interface and the point at which the scattered field will be computed.

    std::vector<std::vector<double>> mesh_pic(N_x + 1, std::vector<double>(N_y + 1, 0.0));
    for (int ii = 0; ii <= N_x; ++ii)
    {
        for (int jj = 0; jj <= N_y; ++jj)
        {
            if (C_Ex[ii][jj] == 0)
                mesh_pic[ii][jj] = N_x / 2; // To get vertical scale in plot OK when plotted
            else if (ii == L)
                mesh_pic[ii][jj] = N_x / 4;
            else if (ii == point1_x && jj == point1_y)
                mesh_pic[ii][jj] = N_x / 2;
        }
    }
    /*
    mesh((1:1:N_y+1),(1:1:N_x+1),mesh_pic);
    axis image;
    title('Simulation region')
    disp('Press any key to continue')
    pause;
    */

    // First time step - Initialize values for H_z, E_x and E_y
    std::vector<std::vector<double>> H_z_nmin1(N_x, std::vector<double>(N_y, 0.0));
    std::vector<std::vector<double>> E_x_nmin1(N_x + 1, std::vector<double>(N_y + 1, 0.0));
    std::vector<std::vector<double>> E_y_nmin1(N_x + 1, std::vector<double>(N_y + 1, 0.0));

    // Pre-allocation
    std::vector<std::vector<double>> H_z_n(N_x, std::vector<double>(N_y, 0.0));
    std::vector<std::vector<double>> E_x_n(N_x + 1, std::vector<double>(N_y + 1, 0.0));
    std::vector<std::vector<double>> E_y_n(N_x + 1, std::vector<double>(N_y + 1, 0.0));

    // Time loop
    for (int m = 1; m < M; ++m)
    {

        // ---------------------------- H field update -----------------------------------------------------------------
        // Campo magnético H_z
        for (int i = 0; i < N_x; ++i)
        {
            for (int j = 0; j < N_y; ++j)
            {
                H_z_n[i][j] += D_Hz * (E_x_nmin1[i][j + 1] - E_x_nmin1[i][j] + E_y_nmin1[i][j] - E_y_nmin1[i + 1][j]);
            }
        }
        // Drive a test line source - used to check basic operation
        // H_z_n(N_x/2,N_y/2) = gaussder((m-1)*delta_t,m_offset,sigma);

        // Special update on scat/tot field boundary
        std::vector<double> E_y_nmin1_inc_front(N_y + 1, Peak * gaussder_norm((m - 1) * delta_t - (L - 1) * delta_s / c, m_offset, sigma));
        // The H_z field is the total field.
        for (int j = L1 - 1; j < N_y - L1; ++j) // talvez deva ser L-1 o indice e iniciar em L1 - 1
            H_z_n[L][j] += H_z_n[L][j] + D_Hz * (E_x_nmin1[L][j + 1] - E_x_nmin1[L][j] + E_y_nmin1[L][j] + E_y_nmin1_inc_front[j] - E_y_nmin1[L + 1][j]);

        // Special update on additional new scat/tot field boundary (only needed for Ey)
        // on the right hand side of Fig. 3.1, at N_x - L1
        // Note now the the "far" side of the ABC is now the SCATTERED, not TOTAL,
        // field, so the role of the E_y fields swops around.
        std::vector<double> E_y_nmin1_inc_back(N_y + 1, Peak * gaussder_norm((m - 1) * delta_t - (N_x - L1 + 1 - 1) * delta_s / c, m_offset, sigma));

        // E_y_nmin1_inc can be overwritten since it is not used again.   Again,
        // the H_z field is the total field.
        for (int j = L1 - 1; j < N_y - L1; ++j)
        {
            H_z_n[N_x - L1][j] = H_z_nmin1[N_x - L1][j] + D_Hz * (E_x_nmin1[N_x - L1][j + 1] - E_x_nmin1[N_x - L1][j] + E_y_nmin1[N_x - L1][j] - E_y_nmin1_inc_back[j] - E_y_nmin1[N_x - L1 + 1][j]);
        }
        // Special update on additional new scat/tot field boundary (only needed for Ey)
        // on the upper side of Fig. 3.1, at N_y - L1. This plays the same role as
        // the ABC at N_x-L1.
        std::vector<double> E_x_nmin1_inc_top(N_y + 1, 0.0); // For this specific polarization (Ey,Hz) the incident x-field is zero

        // The H_z field that follows is the total field. The E_x field above the interface is the scattered field, that
        // below, the total field.
        // Note: since the incident field is zero, the following code stub has not
        // thus been properly tested.
        for (int i = L - 1; i < N_x - L1; ++i)
            H_z_n[i][N_y - L1] += D_Hz * (E_x_nmin1[i][N_y - L1 + 1] + E_x_nmin1_inc_top[i] - E_x_nmin1[i][N_y - L1] + E_y_nmin1[i][N_y - L1] - E_y_nmin1[i + 1][N_y - L1]);

        // Ditto at scat/tot field boundary (only needed for Ey)
        // on the lower side of Fig. 3.1, at L1
        std::vector<double> E_x_nmin1_inc_bottom(N_x + 1, 0.0); // Again, for this specific polarization (Ey,Hz) the incident x-field is zero
        // Note: same caution as above. Also note that role of E_x fields swops
        // around - E_x field above interface is now total field.

        for (int i = L - 1; i < N_x - L1; ++i)
            H_z_n[i][L1] = H_z_nmin1[i][L1] + D_Hz * (E_x_nmin1[i][L1 + 1] - E_x_nmin1[i][L1] - E_x_nmin1_inc_bottom[i] + E_y_nmin1[i][L1] - E_y_nmin1[i + 1][L1]);

        // Very special update - 4 single points on corners of scat/tot field boundary
        // The H_z field is the total field.
        // bottom front corner
        H_z_n[L][L1] = H_z_nmin1[L][L1] + D_Hz * (E_x_nmin1[L][L1 + 1] - E_x_nmin1_inc_bottom[L] - E_x_nmin1[L][L1] + E_y_nmin1[L][L1] + E_y_nmin1_inc_front[L1] - E_y_nmin1[L + 1][L1]);
        // top front corner
        H_z_n[L][N_y - L1] = H_z_nmin1[L][N_y - L1] + D_Hz * (E_x_nmin1[L][N_y - L1 + 1] + E_x_nmin1_inc_top[L] - E_x_nmin1[L][N_y - L1] + E_y_nmin1[L][N_y - L1] + E_y_nmin1_inc_front[N_y - L1] - E_y_nmin1[L + 1][N_y - L1]);

        // At present, following two stubs are commented out - they inject a spurious signal at these corners. It is not clear why.

        /*
         // back bottom corner
         H_z_n(N_x-L1,L1) = H_z_nmin1(N_x-L1,L1) ...
             + D_Hz(N_x-L1,L1).*(  E_x_nmin1(N_x-L1,L1+1) - E_x_nmin1_inc_bottom(N_x-L1)  - E_x_nmin1(N_x-L1,L1) ...
                               + E_y_nmin1(N_x-L1,L1) - E_y_nmin1_inc_back(L1) - E_y_nmin1(N_x-L+1,L1)) ;
         % back top corner
         H_z_n(N_x-L1,N_y-L1) = H_z_nmin1(N_x-L1,N_y-L1) ...
             + D_Hz(N_x-L1,N_y-L1).*(  E_x_nmin1(N_x-L1,N_y-L1+1) + E_x_nmin1_inc_top(N_x-L1)  - E_x_nmin1(N_x-L1,N_y-L1) ...
                               + E_y_nmin1(N_x-L1,N_y-L1) - E_y_nmin1_inc_back(N_y-L1) - E_y_nmin1(N_x-L+1,N_y-L1)) ;
      */

        // ---------------------------- End H field update -----------------------------------------------------------------
        // ---------------------------- E field update -----------------------------------------------------------------
        // Update E fields: (note that latest H fields must be used!)
        // Campo elétrico E_x
        for (int i = 0; i < N_x; ++i)
        {
            for (int j = 1; j < N_y; ++j)
            {
                E_x_n[i][j] += delta_t / (eps_0 * delta_s) * (H_z_n[i][j] - H_z_n[i][j - 1]);
            }
        }

        // Campo elétrico E_y
        for (int i = 1; i < N_x; ++i)
        {
            for (int j = 0; j < N_y; ++j)
            {
                E_y[i][j] -= delta_t / (eps_0 * delta_s) * (H_z[i][j] - H_z[i - 1][j]);
            }
        }

        // Fonte Gaussiana em E_y
        double inc = Peak * gaussder_norm((m - 1) * delta_t - (L - 1) * delta_s / c, m_offset, sigma);
        E_y[L][point1_y] += inc;

        // Salvar valor em ponto específico
        E_y_point1[m] = E_y[point1_x][point1_y];
    }

    std::cout << "Simulação finalizada. Valores de E_y armazenados." << std::endl;

    // Exportar resultados para CSV
    std::ofstream file(out_dir + "/ey_point1.csv");
    file << "tempo_ns,Ey_Vpm\n";
    for (int m = 0; m < M; ++m)
    {
        double time_ns = m * delta_t * 1e9;
        file << time_ns << "," << E_y_point1[m] << "\n";
    }
    file.close();
    std::cout << "Arquivo salvo em: " << out_dir << "/ey_point1.csv" << std::endl;

    return 0;
}