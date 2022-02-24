#include "mesh_tools.hpp"

#include <Galerkin.hpp>

namespace Elasticity
{

#if ELASTICITY_MAX_ELEMENT_ORDER > 2
template <>
void fill_internal_coordinates(MeshType<3> &mesh)
{
    using namespace Galerkin::Transforms;
    for (size_t eli = 0; eli < mesh.num_elements(); ++eli)
    {
        const auto &el_info = mesh.element(eli);
        const auto &nodes = el_info.control_nodes;
        auto transform = TriangleTransform<double>(
            mesh.coord(nodes[0]), mesh.coord(nodes[1]), mesh.coord(nodes[2]));
        mesh.coord_ref(el_info.internal_nodes[0]) = transform(std::array{-0.5, -0.5});
    }
}
#endif

#if ELASTICITY_MAX_ELEMENT_ORDER > 3
template <>
void fill_internal_coordinates(MeshType<4> &mesh)
{
    using namespace Galerkin::Transforms;

    for (size_t eli = 0; eli < mesh.num_elements(); ++eli)
    {
        const auto &el_info = mesh.element(eli);
        const auto &nodes = el_info.control_nodes;
        auto transform = TriangleTransform<double>(
            mesh.coord(nodes[0]), mesh.coord(nodes[1]), mesh.coord(nodes[2]));
        mesh.coord_ref(el_info.internal_nodes[0]) = transform(std::array{-0.5, -0.5});
        mesh.coord_ref(el_info.internal_nodes[1]) = transform(std::array{-0.5, 0.0});
        mesh.coord_ref(el_info.internal_nodes[2]) = transform(std::array{0.0, -0.5});
    }
}
#endif

template <>
void fill_internal_coordinates(MeshVariant &mv)
{
    std::visit([](const auto &mesh) { fill_internal_coordinates(mesh); }, mv);
}

} // namespace Elasticity
