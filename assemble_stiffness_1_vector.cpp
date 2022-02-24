#include "elasticity.hpp"

namespace Elasticity
{

template StiffnessType
assemble_stiffness<Mesh1>(const Mesh1 &, const Eigen::VectorXd &, const Eigen::VectorXd &);

} // namespace Elasticity

