#ifndef TEEIG_ERR_H
#define TEEIG_ERR_H

#include <vector>

std::vector<double> TEeig_err(double a, double b,
                              const std::vector<double>& TEeigvalues,
                              int N, int num_zero_eig);

#endif
