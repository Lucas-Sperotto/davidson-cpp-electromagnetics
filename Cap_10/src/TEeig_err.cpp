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
    std::vector<double> TEexact(k_c.begin(), k_c.begin() + N);

    std::vector<double> rel_err(N);
    for (int i = 0; i < N; ++i) {
        rel_err[i] = std::abs((TEexact[i] - TEeigvalues[num_zero_eig + i]) / TEexact[i]);
    }

    return rel_err;
}
