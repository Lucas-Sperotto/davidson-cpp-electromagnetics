#ifndef CAP_11_EIG_ERR3D_HPP
#define CAP_11_EIG_ERR3D_HPP

#include <vector>

std::vector<double> eig_err3D(double a,
                              double b,
                              double d,
                              const std::vector<double> &eigvalues,
                              int N,
                              int num_zero_eig);

#endif
