#if ELASTICITY_MAX_ELEMENT_ORDER > 3

#include "integrate_traction_force.hpp"

namespace Elasticity
{

template Eigen::Matrix2Xd
integrate_traction_force(const Mesh4 &, size_t, const Eigen::Matrix2Xd &);

} // namespace Elasticity

#endif
