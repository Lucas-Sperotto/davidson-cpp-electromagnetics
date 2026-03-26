#include <vector>
#include <iostream>
#include <cmath>
#include <algorithm> // std::count
#include <tuple>
#include <fstream>
#include <iomanip>
#include <filesystem>
#include <string>
#include <stdexcept>

#include "globals.h" //OK
#include "trimesh.h"
#include "s_nodal.h" //OK
#include "prescr_nodes_mstrip.h"
#include "free_nodes_mstrip.h"

#include <Eigen/Sparse>
#include <Eigen/Dense>

namespace
{
struct Static2DConfig
{
    double a = 2.5;
    double b = 1.25;
    double h = 0.50;
    double w = 0.50;
    double eps_r_sub = 2.55;
    double voltage = 1.0;
    int x_mesh = 20;
    int y_mesh = 20;
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
        << "Uso: ./build/Static2D [opcoes]\n"
        << "  --a VALOR           largura total do problema original (padrao 2.5)\n"
        << "  --b VALOR           altura da caixa (padrao 1.25)\n"
        << "  --h VALOR           altura do substrato (padrao 0.5)\n"
        << "  --w VALOR           largura do condutor central (padrao 0.5)\n"
        << "  --eps-r-sub VALOR   permissividade relativa do substrato (padrao 2.55)\n"
        << "  --voltage VALOR     tensao prescrita no strip (padrao 1.0)\n"
        << "  --x-mesh N          divisao da malha em x (padrao 20)\n"
        << "  --y-mesh N          divisao da malha em y (padrao 20)\n"
        << "  --help              mostra esta ajuda\n";
}
} // namespace

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

int main(int argc, char **argv)
{
    const std::filesystem::path out_dir = PROJECT_OUT_DIR;
    std::filesystem::create_directories(out_dir);

    Static2DConfig config;
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
            if (arg == "--a")
                config.a = std::stod(require_value(i, argc, argv));
            else if (arg == "--b")
                config.b = std::stod(require_value(i, argc, argv));
            else if (arg == "--h")
                config.h = std::stod(require_value(i, argc, argv));
            else if (arg == "--w")
                config.w = std::stod(require_value(i, argc, argv));
            else if (arg == "--eps-r-sub")
                config.eps_r_sub = std::stod(require_value(i, argc, argv));
            else if (arg == "--voltage")
                config.voltage = std::stod(require_value(i, argc, argv));
            else if (arg == "--x-mesh")
                config.x_mesh = std::stoi(require_value(i, argc, argv));
            else if (arg == "--y-mesh")
                config.y_mesh = std::stoi(require_value(i, argc, argv));
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

    if (config.a <= 0.0 || config.b <= 0.0 || config.h <= 0.0 || config.w <= 0.0 ||
        config.eps_r_sub <= 0.0 || config.voltage == 0.0 || config.x_mesh <= 0 || config.y_mesh <= 0)
    {
        std::cerr << "Erro: parametros geometricos, material e malha devem ser positivos; voltage deve ser nao nula.\n";
        return 1;
    }
    if (config.h > config.b)
    {
        std::cerr << "Erro: h nao pode exceder b.\n";
        return 1;
    }

    const double eps_0 = 8.854e-12;
    const double a = config.a;
    const double b = config.b;
    const double h = config.h;
    const double w = config.w;
    const double eps_r_sub = config.eps_r_sub;
    const double V = config.voltage;
    const int x_mesh = config.x_mesh;
    const int y_mesh = config.y_mesh;

    std::cout << "Parametros do Static2D:\n"
              << "  a=" << a << ", b=" << b << ", h=" << h << ", w=" << w << "\n"
              << "  eps_r_sub=" << eps_r_sub << ", V=" << V
              << ", x_mesh=" << x_mesh << ", y_mesh=" << y_mesh << "\n";

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
    Eigen::VectorXd phi_tot;
    if (cg.info() == Eigen::Success)
    {
        phi_tot = cg.solve(b_vec);
        std::cout << "[CG] iterações=" << cg.iterations()
                  << " | erro=" << cg.error() << "\n";
    }
    else
    {
        Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower | Eigen::Upper, Eigen::DiagonalPreconditioner<double>> cg_alt;

        cg_alt.setMaxIterations(max_iter);
        cg_alt.setTolerance(tol);
        cg_alt.compute(A);
        phi_tot = cg_alt.solve(b_vec);
        std::cout << "Fatoracao com IncompleteCholesky falhou; usando pre-condicionador diagonal.\n";
        std::cout << "[CG alt] iterações=" << cg_alt.iterations()
                  << " | erro=" << cg_alt.error() << "\n";
    }

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
        std::ofstream f((out_dir / "field_map.csv").string());
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
        std::ofstream g((out_dir / "strip_nodes.csv").string());
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
        std::ofstream p((out_dir / "geom_params.csv").string());
        p << "a_half,b,h,w,Nx,Ny,eps_r_sub,voltage,capacitance_per_m\n";
        p.setf(std::ios::scientific);
        p << std::setprecision(10);
        p << Lx << "," << Ly << "," << h << "," << w << "," << Nx << "," << Ny << ","
          << eps_r_sub << "," << V << "," << C << "\n";
        p.close();

        std::ofstream summary((out_dir / "static2d_summary.csv").string());
        summary << "key,value\n";
        summary << "a," << a << "\n";
        summary << "b," << b << "\n";
        summary << "h," << h << "\n";
        summary << "w," << w << "\n";
        summary << "eps_r_sub," << eps_r_sub << "\n";
        summary << "voltage," << V << "\n";
        summary << "x_mesh," << x_mesh << "\n";
        summary << "y_mesh," << y_mesh << "\n";
        summary << "num_nodes," << NUM_NODES << "\n";
        summary << "num_elements," << NUM_ELEMENTS << "\n";
        summary << "capacitance_per_m," << C << "\n";
        summary.close();

        std::cout << "Arquivos salvos: field_map.csv, strip_nodes.csv, geom_params.csv, static2d_summary.csv\n";
    }
    // ==================== FIM do bloco ====================
    return 0;
}
