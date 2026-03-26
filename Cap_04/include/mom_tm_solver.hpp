#ifndef CAP_04_MOM_TM_SOLVER_HPP
#define CAP_04_MOM_TM_SOLVER_HPP

#include <complex>
#include <vector>

struct MomTmResult
{
    std::vector<std::complex<double>> current;
    std::vector<double> phi_c;
    double strip_width = 0.0;
    std::vector<double> x_c;
    std::vector<double> y_c;
    double condition_number = 0.0;
};

MomTmResult mom_tm_solver(double k, int N, double a, double E_0, double phi_inc,
                          bool quadrature, double eta, bool toeplitz_flag);

#endif // CAP_04_MOM_TM_SOLVER_HPP
