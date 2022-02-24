#include "elasticity.hpp"

namespace Elasticity
{

namespace detail
{
template <class Mesh>
void add_traction_force(
    const Mesh &mesh, Eigen::Matrix2Xd &nodal_forcing, size_t bindex,
    const Eigen::Matrix2Xd &traction)
{
    const auto &brep = mesh.boundary(bindex);
    size_t bnode = 0;
    for (size_t node : brep.nodes)
    {
        nodal_forcing.col(node) += traction.col(bnode++);
    }
}

} // namespace detail

void add_traction_force(
    const MeshVariant &mv, Eigen::Matrix2Xd &nodal_forcing, size_t bindex,
    const Eigen::Matrix2Xd &traction)
{
    std::visit(
        [&, bindex](const auto &mesh)
        { detail::add_traction_force(mesh, nodal_forcing, bindex, traction); },
        mv);
}

} // namespace Elasticity

