#ifndef RENUMBER_DOF_LTQN_H
#define RENUMBER_DOF_LTQN_H

#include <vector>

void renumber_dof_LTQN(const std::vector<int>& dof_free_flag,
                       std::vector<int>& dof_e1,
                       std::vector<int>& dof_e2,
                       std::vector<int>& dof_f1,
                       std::vector<int>& dof_f2);

#endif
