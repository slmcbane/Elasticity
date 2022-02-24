#include "elasticity.hpp"

namespace Elasticity
{

template <class Mesh>
void impose_dirichlet_condition(
    const Mesh &mesh, StiffnessType &K, Eigen::Matrix2Xd &forcing, size_t which,
    const Eigen::Matrix2Xd &value, double scale)
{
    std::vector<size_t> adjacent;
    adjacent.reserve(2 * max_node_adjacencies);
    const auto &boundary = mesh.boundary(which);

    if (static_cast<size_t>(value.size()) != 2 * boundary.nodes.size())
    {
        throw OutOfBoundsIndex("Vector given for Dirichlet condition has wrong number of values");
    }

    size_t i = 0;
    auto flat_forcing = forcing.reshaped();
    for (auto n : boundary.nodes)
    {
        adjacent.clear();
        for (auto n2 : mesh.adjacent_nodes(n))
        {
            adjacent.push_back(2 * n2);
            adjacent.push_back(2 * n2 + 1);
        }
        K.eliminate_dof(2 * n, value(0, i), scale, flat_forcing, adjacent);
        K.eliminate_dof(2 * n + 1, value(1, i), scale, flat_forcing, adjacent);
        i += 1;
    }
}

} // namespace Elasticity
