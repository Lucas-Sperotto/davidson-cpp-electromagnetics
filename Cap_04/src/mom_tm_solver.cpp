#include "../include/mom_tm_solver.hpp"

#include <Eigen/Dense>

#include <cmath>
#include <complex>
#include <limits>
#include <numbers>
#include <stdexcept>
#include <vector>

namespace
{
std::complex<double> hankel2_0(double x)
{
    return {std::cyl_bessel_j(0.0, x), -std::cyl_neumann(0.0, x)};
}

double condition_number(const Eigen::MatrixXcd &matrix)
{
    Eigen::JacobiSVD<Eigen::MatrixXcd> svd(matrix);
    const auto singular_values = svd.singularValues();
    if (singular_values.size() == 0)
        return 0.0;

    const double sigma_max = singular_values.maxCoeff();
    const double sigma_min = singular_values.minCoeff();
    if (sigma_min <= 0.0)
        return std::numeric_limits<double>::infinity();
    return sigma_max / sigma_min;
}

std::complex<double> efie_kernel(double t, double a, double omega_n, double phi_c,
                                 double delta_phi, double w, double x_c_m,
                                 double y_c_m, double k)
{
    const std::vector<double> tangent = {std::cos(omega_n), std::sin(omega_n)};
    const double phi_1 = phi_c - delta_phi / 2.0;
    const double a1 = std::sqrt(a * a + (w / 2.0) * (w / 2.0));
    const double x = a1 * std::cos(phi_1) + tangent[0] * t;
    const double y = a1 * std::sin(phi_1) + tangent[1] * t;
    const double r_m = std::sqrt((x_c_m - x) * (x_c_m - x) + (y_c_m - y) * (y_c_m - y));
    return hankel2_0(k * r_m);
}
} // namespace

MomTmResult mom_tm_solver(double k, int N, double a, double E_0, double phi_inc,
                          bool quadrature, double eta, bool toeplitz_flag)
{
    // Translation of MoM_TM_solver.m.
    if (N <= 0)
        throw std::invalid_argument("N deve ser positivo.");
    if (a <= 0.0 || k <= 0.0 || eta <= 0.0)
        throw std::invalid_argument("a, k e eta devem ser positivos.");

    MomTmResult result;
    result.current.assign(N, {0.0, 0.0});
    result.phi_c.assign(N, 0.0);
    result.x_c.assign(N, 0.0);
    result.y_c.assign(N, 0.0);

    Eigen::VectorXcd V_vec = Eigen::VectorXcd::Zero(N);
    Eigen::MatrixXcd Z_mat = Eigen::MatrixXcd::Zero(N, N);

    const double delta_phi = 2.0 * std::numbers::pi / static_cast<double>(N);
    std::vector<double> omega(N, 0.0);

    // Pre-compute segment midpoint angles and segment orientation angles.
    for (int mm = 0; mm < N; ++mm)
    {
        result.phi_c[mm] = mm * delta_phi;
        omega[mm] = std::numbers::pi / 2.0 + mm * delta_phi;
    }

    // Width of strips. The polygon circumscribes the circular cylinder,
    // matching the assumption used in the MATLAB code.
    result.strip_width = 2.0 * a * std::tan(delta_phi / 2.0);
    for (int mm = 0; mm < N; ++mm)
    {
        result.x_c[mm] = a * std::cos(result.phi_c[mm]);
        result.y_c[mm] = a * std::sin(result.phi_c[mm]);
    }

    const int m_stop = toeplitz_flag ? 1 : N;

    for (int mm = 0; mm < m_stop; ++mm)
    {
        for (int nn = 0; nn < N; ++nn)
        {
            if (mm == nn)
            {
                constexpr double gamma = 1.781072418;
                const std::complex<double> imag_unit(0.0, 1.0);
                Z_mat(mm, nn) = k * eta * result.strip_width / 4.0 *
                                (1.0 - imag_unit * 2.0 / std::numbers::pi *
                                           (std::log(gamma * k * result.strip_width / 4.0) - 1.0));
            }
            else if (!quadrature)
            {
                const double r_mn =
                    std::sqrt((result.x_c[mm] - result.x_c[nn]) * (result.x_c[mm] - result.x_c[nn]) +
                              (result.y_c[mm] - result.y_c[nn]) * (result.y_c[mm] - result.y_c[nn]));
                // Approximate eq. 2.30 using a single midpoint evaluation.
                Z_mat(mm, nn) = k * eta / 4.0 * result.strip_width * hankel2_0(k * r_mn);
            }
            else
            {
                // Evaluate eq. 2.8 with a simple low-order quadrature.
                constexpr int n_quad = 10;
                const double delta_t = result.strip_width / static_cast<double>(n_quad);
                std::complex<double> z_entry = 0.0;
                for (int ii = 0; ii < n_quad; ++ii)
                {
                    const double t = (ii + 0.5) * delta_t;
                    z_entry += efie_kernel(t, a, omega[nn], result.phi_c[nn], delta_phi,
                                           result.strip_width, result.x_c[mm], result.y_c[mm], k) *
                               delta_t;
                }
                Z_mat(mm, nn) = k * eta / 4.0 * z_entry;
            }
        }
    }

    if (toeplitz_flag)
    {
        const Eigen::VectorXcd first_row = Z_mat.row(0);
        for (int mm = 0; mm < N; ++mm)
        {
            for (int nn = 0; nn < N; ++nn)
            {
                Z_mat(mm, nn) = first_row(std::abs(nn - mm));
            }
        }
    }

    // Incident field sampled at segment centres.
    for (int mm = 0; mm < N; ++mm)
    {
        V_vec(mm) =
            E_0 * std::exp(std::complex<double>(0.0, 1.0) *
                           k * (result.x_c[mm] * std::cos(phi_inc) + result.y_c[mm] * std::sin(phi_inc)));
    }

    const Eigen::VectorXcd I_vec = Z_mat.partialPivLu().solve(V_vec);
    result.condition_number = condition_number(Z_mat);
    for (int mm = 0; mm < N; ++mm)
        result.current[mm] = I_vec(mm);

    return result;
}
