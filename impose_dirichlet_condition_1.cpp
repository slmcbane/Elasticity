#include "impose_dirichlet_condition.hpp"

namespace Elasticity
{

template void impose_dirichlet_condition(
    const Mesh1 &, StiffnessType &, Eigen::Matrix2Xd &, size_t, const Eigen::Matrix2Xd &, double);

} // namespace Elasticity
