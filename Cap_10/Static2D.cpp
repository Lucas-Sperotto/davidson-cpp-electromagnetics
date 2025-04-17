#include <vector>
#include <iostream>
#include <cmath>

int NUM_NODES, NUM_ELEMENTS, NUM_DOFS;
std::vector<std::vector<int>> ELEMENTS;
std::vector<std::vector<double>> NODE_COORD;

void Static2D() {
    const double eps_0 = 8.854e-12;
    const double mu_0 = 4.0 * M_PI * 1e-7;

    double a = 2.5;
    double b = a / 2.0;
    double h = 0.5;
    double w = h;
    double eps_r_sub = 2.55;
    double V = 1.0;

    int x_mesh = 20, y_mesh = 20;
    std::vector<std::vector<double>> x_nodes, y_nodes;

    // Gera malha e preenche ELEMENTS e NODE_COORD
    trimesh(a / 2.0, b, x_mesh, y_mesh, x_nodes, y_nodes);

    // Centro de cada elemento para uso posterior
    std::vector<double> x_c(NUM_ELEMENTS), y_c(NUM_ELEMENTS);

    for (int i = 0; i < NUM_NODES; ++i)
        std::cout << "Node " << i << ": (" << NODE_COORD[i][0] << ", " << NODE_COORD[i][1] << ")\n";

    for (int e = 0; e < NUM_ELEMENTS; ++e) {
        double xc = 0, yc = 0;
        for (int j = 0; j < 3; ++j) {
            xc += NODE_COORD[ELEMENTS[e][j]][0];
            yc += NODE_COORD[ELEMENTS[e][j]][1];
        }
        x_c[e] = xc / 3.0;
        y_c[e] = yc / 3.0;
        std::cout << "Elem " << e << ": (" << x_c[e] << ", " << y_c[e] << ")\n";
    }

    // Definição dos materiais
    std::vector<double> eps_r(NUM_ELEMENTS, 1.0);
    for (int e = 0; e < NUM_ELEMENTS; ++e) {
        if (y_c[e] < h) eps_r[e] = eps_r_sub;
    }

    std::cout << "Malha e materiais definidos. Pressione ENTER para continuar.\n";
    std::cin.get();


    #include <algorithm>  // std::count

    // Externas
    extern std::vector<std::vector<double>> NODE_COORD;
    extern int NUM_NODES;
    
    // Funções auxiliares (pré-definidas)
    void free_nodes_mstrip(double a, double b, double h, double w,
                           std::vector<int>& node_free_flag, int& num_free_nodes);
    void prescr_nodes_mstrip(double a, double b, double h, double w,
                             std::vector<int>& node_prenz_flag, int& num_prenz_nodes);
    
    ...
    
        // Marca nós livres
        std::vector<int> node_free_flag;
        int num_free_nodes = 0;
        free_nodes_mstrip(a, b, h, w, node_free_flag, num_free_nodes);
    
        // Marca nós prescritos (complemento)
        std::vector<int> node_pre_flag(NUM_NODES);
        for (int i = 0; i < NUM_NODES; ++i)
            node_pre_flag[i] = !node_free_flag[i];
    
        // Marca nós com valor prescrito (condutor central)
        std::vector<int> node_prenz_flag;
        int num_prenz_nodes = 0;
        prescr_nodes_mstrip(a, b, h, w, node_prenz_flag, num_prenz_nodes);
    
        // Inicializa vetor de potenciais prescritos
        std::vector<double> phi_pre(NUM_NODES, 0.0);
        for (int i = 0; i < NUM_NODES; ++i) {
            if (node_prenz_flag[i])
                phi_pre[i] = V;
        }
    
#include <Eigen/Sparse>
#include <Eigen/Dense>

// s_nodal deve estar implementada: retorna matriz local 3x3
std::vector<std::vector<double>> s_nodal(double x1, double y1,
                                         double x2, double y2,
                                         double x3, double y3);

...

    // Inicialização do sistema
    Eigen::SparseMatrix<double> S_mat(NUM_NODES, NUM_NODES);
    Eigen::VectorXd b_vec = Eigen::VectorXd::Zero(NUM_NODES);

    std::vector<Eigen::Triplet<double>> tripletList;

    for (int ielem = 0; ielem < NUM_ELEMENTS; ++ielem) {
        const auto& tri = ELEMENTS[ielem];

        double x1 = NODE_COORD[tri[0]][0], y1 = NODE_COORD[tri[0]][1];
        double x2 = NODE_COORD[tri[1]][0], y2 = NODE_COORD[tri[1]][1];
        double x3 = NODE_COORD[tri[2]][0], y3 = NODE_COORD[tri[2]][1];

        auto S_elem = s_nodal(x1, y1, x2, y2, x3, y3);
        for (auto& row : S_elem)
            for (auto& val : row)
                val *= eps_r[ielem]; // material

        for (int j = 0; j < 3; ++j) {
            int jj = tri[j];
            for (int k = 0; k < 3; ++k) {
                int kk = tri[k];

                if (node_free_flag[jj] && node_free_flag[kk]) {
                    tripletList.emplace_back(jj, kk, S_elem[j][k]);
                }
                else if (node_free_flag[jj] && node_pre_flag[kk]) {
                    b_vec(jj) -= S_elem[j][k] * phi_pre[kk];

                    // Para evitar singularidades, força valor de φ prescrito na diagonal
                    tripletList.emplace_back(kk, kk, 1.0);
                    b_vec(kk) = phi_pre[kk];
                }
            }
        }
    }

    S_mat.setFromTriplets(tripletList.begin(), tripletList.end());


    #include <Eigen/Sparse>
    #include <Eigen/Dense>
    
    // s_nodal deve estar implementada: retorna matriz local 3x3
    std::vector<std::vector<double>> s_nodal(double x1, double y1,
                                             double x2, double y2,
                                             double x3, double y3);
    
    ...
    
        // Inicialização do sistema
        Eigen::SparseMatrix<double> S_mat(NUM_NODES, NUM_NODES);
        Eigen::VectorXd b_vec = Eigen::VectorXd::Zero(NUM_NODES);
    
        std::vector<Eigen::Triplet<double>> tripletList;
    
        for (int ielem = 0; ielem < NUM_ELEMENTS; ++ielem) {
            const auto& tri = ELEMENTS[ielem];
    
            double x1 = NODE_COORD[tri[0]][0], y1 = NODE_COORD[tri[0]][1];
            double x2 = NODE_COORD[tri[1]][0], y2 = NODE_COORD[tri[1]][1];
            double x3 = NODE_COORD[tri[2]][0], y3 = NODE_COORD[tri[2]][1];
    
            auto S_elem = s_nodal(x1, y1, x2, y2, x3, y3);
            for (auto& row : S_elem)
                for (auto& val : row)
                    val *= eps_r[ielem]; // material
    
            for (int j = 0; j < 3; ++j) {
                int jj = tri[j];
                for (int k = 0; k < 3; ++k) {
                    int kk = tri[k];
    
                    if (node_free_flag[jj] && node_free_flag[kk]) {
                        tripletList.emplace_back(jj, kk, S_elem[j][k]);
                    }
                    else if (node_free_flag[jj] && node_pre_flag[kk]) {
                        b_vec(jj) -= S_elem[j][k] * phi_pre[kk];
    
                        // Para evitar singularidades, força valor de φ prescrito na diagonal
                        tripletList.emplace_back(kk, kk, 1.0);
                        b_vec(kk) = phi_pre[kk];
                    }
                }
            }
        }
    
        S_mat.setFromTriplets(tripletList.begin(), tripletList.end());
    



        double compute_static_capacitance(const Eigen::VectorXd& phi_tot,
            const std::vector<double>& eps_r,
            double eps_0, double V) {
double E_tot = 0.0;

for (int ielem = 0; ielem < NUM_ELEMENTS; ++ielem) {
const auto& tri = ELEMENTS[ielem];

double x1 = NODE_COORD[tri[0]][0], y1 = NODE_COORD[tri[0]][1];
double x2 = NODE_COORD[tri[1]][0], y2 = NODE_COORD[tri[1]][1];
double x3 = NODE_COORD[tri[2]][0], y3 = NODE_COORD[tri[2]][1];

auto S_elem = s_nodal(x1, y1, x2, y2, x3, y3);
for (auto& row : S_elem)
for (auto& val : row)
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
