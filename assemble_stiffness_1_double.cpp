#include "elasticity.hpp"

namespace Elasticity
{

template StiffnessType assemble_stiffness<Mesh1>(const Mesh1 &, double, double);

} // namespace Elasticity

