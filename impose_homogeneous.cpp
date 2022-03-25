#include "elasticity.hpp"

namespace Elasticity
{

template <>
void impose_homogeneous_condition(
    const MeshVariant &mv, StiffnessType &K, size_t which, double scale)
{
    std::visit(
        [&K, which, scale](const auto &mesh)
        { impose_homogeneous_condition(mesh, K, which, scale); },
        mv);
}

} // namespace Elasticity
