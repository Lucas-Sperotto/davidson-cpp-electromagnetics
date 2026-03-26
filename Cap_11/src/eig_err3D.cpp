#include "../include/eig_err3D.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <vector>

// Translation of eig_err3D.m:
// relative error in the first N non-zero eigenvalues for the box test case.
std::vector<double> eig_err3D(double a,
                              double b,
                              double d,
                              const std::vector<double> &eigvalues,
                              int N,
                              int num_zero_eig)
{
    (void)a;
    (void)b;
    (void)d;

    const std::array<double, 8> exact = {
        5.23599, 7.02481, 7.55145, 7.55145,
        8.17887, 8.17887, 8.88577, 8.94726};

    const int count = std::min({N, static_cast<int>(exact.size()),
                                static_cast<int>(eigvalues.size()) - num_zero_eig});
    std::vector<double> rel_err;
    rel_err.reserve(count);
    for (int index = 0; index < count; ++index)
    {
        rel_err.push_back(
            std::abs((exact[index] - eigvalues[num_zero_eig + index]) / exact[index]));
    }
    return rel_err;
}
