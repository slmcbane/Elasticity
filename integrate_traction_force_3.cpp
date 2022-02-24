#if ELASTICITY_MAX_ELEMENT_ORDER > 2

#include "integrate_traction_force.hpp"

namespace Elasticity
{

template Eigen::Matrix2Xd
integrate_traction_force(const Mesh3 &, size_t, const Eigen::Matrix2Xd &);

} // namespace Elasticity

#endif
