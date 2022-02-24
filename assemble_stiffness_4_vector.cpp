#if ELASTICITY_MAX_ELEMENT_ORDER > 3

#include "elasticity.hpp"

namespace Elasticity
{

template StiffnessType
assemble_stiffness<Mesh4>(const Mesh4 &, const Eigen::VectorXd &, const Eigen::VectorXd &);

} // namespace Elasticity

#endif
