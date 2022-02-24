#include "elasticity.hpp"

namespace Elasticity
{

template StiffnessType
assemble_stiffness<Mesh2>(const Mesh2 &, const Eigen::VectorXd &, const Eigen::VectorXd &);

} // namespace Elasticity
