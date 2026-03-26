#include <vector>
#include <cmath>
#include <algorithm>

// Calcula o erro relativo dos N primeiros autovalores não-nulos
std::vector<double> TEeig_err(double a, double b,
                              const std::vector<double>& TEeigvalues,
                              int N, int num_zero_eig) {
    std::vector<double> k_c;
    for (int m = 0; m <= 10; ++m) {
        for (int n = 0; n <= 9; ++n) {
            if (m != 0 || n != 0) {
                double kc = std::sqrt(std::pow(m * M_PI / a, 2) + std::pow(n * M_PI / b, 2));
                k_c.push_back(kc);
            }
        }
    }

    std::sort(k_c.begin(), k_c.end());
    int available_eigs = static_cast<int>(TEeigvalues.size()) - num_zero_eig;
    int usable = std::min<int>(N, std::min<int>(static_cast<int>(k_c.size()), available_eigs));
    if (usable <= 0) {
        return {};
    }

    std::vector<double> TEexact(k_c.begin(), k_c.begin() + usable);

    std::vector<double> rel_err(usable);
    for (int i = 0; i < usable; ++i) {
        rel_err[i] = std::abs((TEexact[i] - TEeigvalues[num_zero_eig + i]) / TEexact[i]);
    }

    return rel_err;
}
