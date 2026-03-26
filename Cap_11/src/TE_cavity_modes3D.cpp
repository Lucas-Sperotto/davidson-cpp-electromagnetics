#include "../include/TE_cavity_modes3D.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <utility>
#include <vector>

// Translation of TE_cavity_modes3D.m:
// first approximately TE_z cavity eigenmodes in a rectangular cavity.
std::pair<std::vector<double>, std::vector<std::array<int, 3>>> TE_cavity_modes3D(
    double a,
    double b,
    double d,
    int max_index)
{
    std::vector<std::pair<double, std::array<int, 3>>> modes;
    for (int mm = 0; mm <= max_index; ++mm)
    {
        for (int nn = 0; nn <= max_index; ++nn)
        {
            for (int pp = 0; pp <= max_index; ++pp)
            {
                if ((mm != 0) + (nn != 0) + (pp != 0) >= 2)
                {
                    const double value = std::sqrt(
                        std::pow(mm * M_PI / a, 2) +
                        std::pow(nn * M_PI / b, 2) +
                        std::pow(pp * M_PI / d, 2));
                    modes.push_back({value, {mm, nn, pp}});
                }
            }
        }
    }

    std::sort(modes.begin(), modes.end(),
              [](const auto &lhs, const auto &rhs)
              { return lhs.first < rhs.first; });

    std::vector<double> eigvalues;
    std::vector<std::array<int, 3>> mode_index;
    eigvalues.reserve(modes.size());
    mode_index.reserve(modes.size());
    for (const auto &entry : modes)
    {
        eigvalues.push_back(entry.first);
        mode_index.push_back(entry.second);
    }
    return {eigvalues, mode_index};
}
