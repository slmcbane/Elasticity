#if ELASTICITY_MAX_ELEMENT_ORDER > 3

#include "impose_dirichlet_condition.hpp"

namespace Elasticity
{

template void impose_dirichlet_condition(
    const Mesh4 &, StiffnessType &, Eigen::Matrix2Xd &, size_t, const Eigen::Matrix2Xd &, double);

} // namespace Elasticity

#endif
