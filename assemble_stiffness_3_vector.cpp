#if ELASTICITY_MAX_ELEMENT_ORDER > 2

#include "elasticity.hpp"

namespace Elasticity
{

template StiffnessType
assemble_stiffness<Mesh3>(const Mesh3 &, const Eigen::VectorXd &, const Eigen::VectorXd &);

} // namespace Elasticity

#endif

