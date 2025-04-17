#ifndef PLOT_FIELD_H
#define PLOT_FIELD_H

#include <vector>

void plot_field(const std::vector<double>& dofs,
                const std::vector<int>& dof_e1,
                const std::vector<double>& XX,
                const std::vector<double>& YY,
                int rows, int cols, int plotnum,
                double eigvalue);

#endif
