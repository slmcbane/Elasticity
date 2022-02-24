#if ELASTICITY_MAX_ELEMENT_ORDER > 2

#include "elasticity.hpp"
#include "integrate_volume_force.hpp"

namespace Elasticity
{

template Eigen::Matrix2Xd integrate_volume_force(const Mesh3 &, const Eigen::Matrix2Xd &);

}

#endif
