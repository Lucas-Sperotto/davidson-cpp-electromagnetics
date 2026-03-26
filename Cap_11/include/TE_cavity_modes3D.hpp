#ifndef CAP_11_TE_CAVITY_MODES3D_HPP
#define CAP_11_TE_CAVITY_MODES3D_HPP

#include <array>
#include <utility>
#include <vector>

std::pair<std::vector<double>, std::vector<std::array<int, 3>>> TE_cavity_modes3D(
    double a,
    double b,
    double d,
    int max_index);

#endif
