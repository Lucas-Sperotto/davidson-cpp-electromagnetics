#include <vector>
#include <iostream>
#include <cmath>
#include <algorithm> // std::count
#include <tuple>
#include <fstream>
#include <iomanip>

#include "globals.h" //OK
#include "trimesh.h"
#include "s_nodal.h" //OK
#include "prescr_nodes_mstrip.h"
#include "free_nodes_mstrip.h"

#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/Dense>

double compute_static_capacitance(const Eigen::VectorXd &phi_tot,
                                  const std::vector<double> &eps_r,
                                  double eps_0, double V)
{
    double E_tot = 0.0;

    for (int ielem = 0; ielem < NUM_ELEMENTS; ++ielem)
    {
        const auto &tri = ELEMENTS[ielem];

        double x1 = NODE_COORD[tri[0]][0], y1 = NODE_COORD[tri[0]][1];
        double x2 = NODE_COORD[tri[1]][0], y2 = NODE_COORD[tri[1]][1];
        double x3 = NODE_COORD[tri[2]][0], y3 = NODE_COORD[tri[2]][1];

        auto S_elem = s_nodal(x1, y1, x2, y2, x3, y3);
        for (auto &row : S_elem)
            for (auto &val : row)
                val *= eps_r[ielem];

        Eigen::Vector3d phi_elem;
        for (int i = 0; i < 3; ++i)
            phi_elem(i) = phi_tot(tri[i]);

        double E_elem = 0.0;
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                E_elem += phi_elem(i) * S_elem[i][j] * phi_elem(j);
        E_tot += 0.5 * E_elem;
    }

    E_tot *= 2.0 * eps_0; // fator 2 pela simetria

    double C = 2.0 * E_tot / (V * V);
    return C;
}

int main()
{
    const double eps_0 = 8.854e-12;
    const double mu_0 = 4.0 * M_PI * 1e-7;
    const double mu_r = 1.0;
    double a = 2.5;
    double b = a / 2.0;
    double h = 0.50;
    double w = h;
    double eps_r_sub = 1.0; // 2.55; // 2.55;//2.55;
    double V = 1.0;
    // std::cout << "Nodes: " << NUM_NODES << "\n";
    int x_mesh = 100, y_mesh = 100;
    std::vector<std::vector<double>> x_nodes, y_nodes;

    // Gera malha e preenche ELEMENTS e NODE_COORD
    trimesh(a / 2.0, b, x_mesh, y_mesh, x_nodes, y_nodes);
    // std::cout << "Nodes: " << NUM_NODES << "\n";
    //  Centro de cada elemento para uso posterior
    std::vector<double> x_c(NUM_ELEMENTS), y_c(NUM_ELEMENTS);

    // for (int i = 0; i < NUM_NODES; ++i)
    //  std::cout << "Node " << i << ": (" << NODE_COORD[i][0] << ", " << NODE_COORD[i][1] << ")\n";

    for (int ielem = 0; ielem < NUM_ELEMENTS; ++ielem)
    {
        double xc = 0, yc = 0;
        for (int j = 0; j < 3; ++j)
        {
            xc += NODE_COORD[ELEMENTS[ielem][j]][0];
            yc += NODE_COORD[ELEMENTS[ielem][j]][1];
        }
        x_c[ielem] = xc / 3.0;
        y_c[ielem] = yc / 3.0;
        // std::cout << "Elem " << ielem << ": (" << x_c[ielem] << ", " << y_c[ielem] << ")\n";
    }

    // Definição dos materiais
    std::vector<double> eps_r(NUM_ELEMENTS, 1.0);
    for (int ielem = 0; ielem < NUM_ELEMENTS; ++ielem)
    {
        if (y_c[ielem] < h)
            eps_r[ielem] = eps_r_sub;
    }

    std::cout << "Malha e materiais definidos.\n";

    // Marca nós livres
    std::vector<int> node_free_flag;
    int num_free_nodes = 0;
    free_nodes_mstrip(a, b, h, w, node_free_flag, num_free_nodes);

    std::cout << "Nós livres definidos.\n";

    // Marca nós prescritos (complemento)
    std::vector<int> node_pre_flag(NUM_NODES);
    for (int i = 0; i < NUM_NODES; ++i)
        node_pre_flag[i] = !node_free_flag[i];

    // Marca nós com valor prescrito (condutor central)
    std::vector<int> node_prenz_flag;
    int num_prenz_nodes = 0;
    prescr_nodes_mstrip(a, b, h, w, node_prenz_flag, num_prenz_nodes);

    std::cout << "Nós prescritos definidos.\n";

    // Inicializa vetor de potenciais prescritos
    std::vector<double> phi_pre(NUM_NODES, 0.0);
    for (int i = 0; i < NUM_NODES; ++i)
    {
        if (node_prenz_flag[i])
            phi_pre[i] = V;
    }

    // Inicialização do sistema (denso)
    Eigen::MatrixXd S_mat = Eigen::MatrixXd::Zero(NUM_NODES, NUM_NODES);
    Eigen::VectorXd b_vec = Eigen::VectorXd::Zero(NUM_NODES);

    std::cout << "Potenciais prescritos definidos.\n";

    std::cout << "iniciando laço dos elementos.\n";
    for (int ielem = 0; ielem < NUM_ELEMENTS; ++ielem)
    {
        const auto &tri = ELEMENTS[ielem];

        double x1 = NODE_COORD[tri[0]][0], y1 = NODE_COORD[tri[0]][1];
        double x2 = NODE_COORD[tri[1]][0], y2 = NODE_COORD[tri[1]][1];
        double x3 = NODE_COORD[tri[2]][0], y3 = NODE_COORD[tri[2]][1];

        auto S_elem = s_nodal(x1, y1, x2, y2, x3, y3);
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                S_elem[i][j] *= eps_r[ielem]; // material

        for (int jnode = 0; jnode < 3; ++jnode)
        {
            int jj = tri[jnode];
            for (int knode = 0; knode < 3; ++knode)
            {
                int kk = tri[knode];

                if (node_free_flag[jj] && node_free_flag[kk]) // livre–livre: entra na matriz
                    S_mat(jj, kk) += S_elem[jnode][knode];
                else if (node_free_flag[jj] && node_pre_flag[kk]) // livre–prescrito: vai para o termo fonte
                    b_vec(jj) -= S_elem[jnode][knode] * phi_pre[kk];
                // casos com jj prescrito são ignorados (não montamos linhas/cols prescritas)
            }
        }
    }

    std::cout << "terminou laço dos elementos..\n";

    // Impor nós prescritos UMA vez (fora da montagem por elemento)
    for (int k = 0; k < NUM_NODES; ++k)
    {
        if (node_pre_flag[k])
        {
            // (opcional, mas seguro em matriz densa)
            S_mat.row(k).setZero();
            S_mat.col(k).setZero();

            S_mat(k, k) = 1.0;
            b_vec(k) = phi_pre[k];
        }
    }

    std::cout << "Imposto condição de contorno.\n";
    std::cout << "Iniciando solver.\n";

    int max_iter = 2000;
    double tol = 1e-10;

    Eigen::SparseMatrix<double> A = S_mat.sparseView(1.0, 1e-16); // (refDensity=1.0, refPrune=prune_tol)
    A.makeCompressed();

    // 3) Conjugate Gradient + pré-condicionador (SPD)
    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower | Eigen::Upper, Eigen::IncompleteCholesky<double>> cg;

    cg.setMaxIterations(max_iter);
    cg.setTolerance(tol);
    cg.compute(A);
    std::cout << "Resolvendo...\n";
    //  Caso a fatoração do pré-condicionador falhe, ainda podemos tentar resolver
    Eigen::VectorXd phi_tot;
    if (cg.info() != Eigen::Success)
    {
        // Tenta sem pré-condicionador (diagonal)
        Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower | Eigen::Upper, Eigen::DiagonalPreconditioner<double>> cg_alt;

        cg_alt.setMaxIterations(max_iter);
        cg_alt.setTolerance(tol);
        cg_alt.compute(A);
        phi_tot = cg_alt.solve(b_vec);
        std::cout << "Resolvendo com pré condicionador...\n";
        std::cout << "[CG] iterações=" << cg_alt.iterations()
                  << " | erro=" << cg_alt.error() << "\n";
    }

    phi_tot = cg.solve(b_vec);
    std::cout << "[CG] iterações=" << cg.iterations()
              << " | erro=" << cg.error() << "\n";

    // Capacitância (sua função existente)
    std::cout << "Solução Calculada.\n";
    double C = compute_static_capacitance(phi_tot, eps_r, eps_0, V);

    std::cout << "Capacitância por unidade de comprimento: C = " << C << " F/m\n";

    // ==================== Pós-solução: campo elétrico e CSV ====================
    {
        // garanta que 'phi_tot' seja o vetor final do solver (sem sombra de variável)
        // se você mantiver os dois caminhos (cg e cg_alt), escolha um para preencher 'phi_tot_final'

        // Parâmetros geométricos e malha tal como usados
        const int Nx = x_mesh;  // número de células em x
        const int Ny = y_mesh;  // número de células em y
        const int nx1 = Nx + 1; // nós em x
        const int ny1 = Ny + 1; // nós em y
        const double Lx = a / 2.0;
        const double Ly = b;
        const double dx = Lx / Nx;
        const double dy = Ly / Ny;

        // (1) Reorganiza phi em grade 2D com a mesma indexação do MATLAB:
        // idx = i + nx1*j  (i = 0..Nx, j = 0..Ny)
        std::vector<double> phi_grid(nx1 * ny1);
        for (int j = 0; j < ny1; ++j)
        {
            for (int i = 0; i < nx1; ++i)
            {
                int idx = i + nx1 * j;
                phi_grid[idx] = phi_tot(idx);
            }
        }

        // (2) Calcula Ex, Ey = -grad(phi)
        std::vector<double> Ex_grid(nx1 * ny1, 0.0);
        std::vector<double> Ey_grid(nx1 * ny1, 0.0);

        auto phi_at = [&](int i, int j) -> double
        {
            return phi_grid[i + nx1 * j];
        };

        for (int j = 0; j < ny1; ++j)
        {
            for (int i = 0; i < nx1; ++i)
            {
                double dphidx, dphidy;

                // dphi/dx
                if (i == 0)
                {
                    dphidx = (phi_at(1, j) - phi_at(0, j)) / dx; // avanço
                }
                else if (i == nx1 - 1)
                {
                    dphidx = (phi_at(nx1 - 1, j) - phi_at(nx1 - 2, j)) / dx; // retrocesso
                }
                else
                {
                    dphidx = (phi_at(i + 1, j) - phi_at(i - 1, j)) / (2.0 * dx); // central
                }

                // dphi/dy
                if (j == 0)
                {
                    dphidy = (phi_at(i, 1) - phi_at(i, 0)) / dy;
                }
                else if (j == ny1 - 1)
                {
                    dphidy = (phi_at(i, ny1 - 1) - phi_at(i, ny1 - 2)) / dy;
                }
                else
                {
                    dphidy = (phi_at(i, j + 1) - phi_at(i, j - 1)) / (2.0 * dy);
                }

                int idx = i + nx1 * j;
                Ex_grid[idx] = -dphidx; // E = -grad(phi)
                Ey_grid[idx] = -dphidy;
            }
        }

        // (3) Salva campo + potencial em CSV
        std::ofstream f("../out/field_map.csv");
        f << "x,y,phi,Ex,Ey\n";
        f.setf(std::ios::scientific);
        f << std::setprecision(10);

        for (int j = 0; j < ny1; ++j)
        {
            for (int i = 0; i < nx1; ++i)
            {
                int idx = i + nx1 * j;
                double x = i * dx;
                double y = j * dy;
                f << x << "," << y << ","
                  << phi_grid[idx] << ","
                  << Ex_grid[idx] << ","
                  << Ey_grid[idx] << "\n";
            }
        }
        f.close();

        // (4) Exporta nós do microstrip (prescritos a V) para desenhar o condutor
        std::ofstream g("../out/strip_nodes.csv");
        g << "x,y\n";
        g.setf(std::ios::scientific);
        g << std::setprecision(10);
        for (int i = 0; i < NUM_NODES; ++i)
        {
            if (node_prenz_flag[i])
            {
                g << NODE_COORD[i][0] << "," << NODE_COORD[i][1] << "\n";
            }
        }
        g.close();

        // (5) Opcional: salva parâmetros para o Python ler (facilita manter coerência)
        std::ofstream p("../out/geom_params.csv");
        p << "a_half,b,h,w,Nx,Ny\n";
        p.setf(std::ios::scientific);
        p << std::setprecision(10);
        p << Lx << "," << Ly << "," << h << "," << w << "," << Nx << "," << Ny << "\n";
        p.close();

        std::cout << "Arquivos salvos: field_map.csv, strip_nodes.csv, geom_params.csv\n";
    }
    // ==================== FIM do bloco ====================
    return 0;
}
