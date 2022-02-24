#include "read_mesh.hpp"
#include "mesh_tools.hpp"
#include <cstdio>

namespace Elasticity
{

MeshVariant read_mesh(const char *name, int order)
{
    assert(order == 1 || order == 2 || order == 3 || order == 4);
    std::optional<MeshVariant> maybe_variant;
    if (order == 1)
    {
        maybe_variant =
            msh::parse_gmsh_to_tetmesh<max_element_adjacencies, max_node_adjacencies, 0, 0>(name);
    }
    else if (order == 2)
    {
        maybe_variant =
            msh::parse_gmsh_to_tetmesh<max_element_adjacencies, max_node_adjacencies, 1, 0>(name);
    }
#if ELASTICITY_MAX_ELEMENT_ORDER > 2
    else if (order == 3)
    {
        maybe_variant =
            msh::parse_gmsh_to_tetmesh<max_element_adjacencies, max_node_adjacencies, 2, 1>(name);
    }
#if ELASTICITY_MAX_ELEMENT_ORDER > 3
    else
    {
        maybe_variant =
            msh::parse_gmsh_to_tetmesh<max_element_adjacencies, max_node_adjacencies, 3, 3>(name);
    }
#endif
#endif

    if (!maybe_variant)
    {
        fprintf(stderr, "Got invalid mesh order %d in read_mesh!\n", order);
        exit(1);
    }

    std::visit(
        [](auto &mesh)
        {
            mesh.renumber_nodes();
            fill_internal_coordinates(mesh);
        },
        maybe_variant.value());

    return *maybe_variant;
}

} // namespace Elasticity
