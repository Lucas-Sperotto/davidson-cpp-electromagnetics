#ifndef LTQN_H
#define LTQN_H

#include <vector>

std::vector<std::vector<double>> LTQN(
    const std::vector<double>& lambda,
    const std::vector<std::vector<double>>& nabla_lambda);

#endif
