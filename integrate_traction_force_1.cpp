#include "integrate_traction_force.hpp"

namespace Elasticity
{

template Eigen::Matrix2Xd
integrate_traction_force(const Mesh1 &, size_t, const Eigen::Matrix2Xd &);

} // namespace Elasticity
