#include "elasticity.hpp"
#include "integrate_volume_force.hpp"

namespace Elasticity
{

template Eigen::Matrix2Xd integrate_volume_force(const Mesh1 &, const Eigen::Matrix2Xd &);

} // namespace Elasticity
