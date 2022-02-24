#include "mesh_tools.hpp"

namespace Elasticity
{

template <>
std::optional<size_t> find_boundary_with_tag(const MeshVariant &mv, std::string_view tag)
{
    return std::visit([tag](const auto &mesh) { return find_boundary_with_tag(mesh, tag); }, mv);
}

} // namespace Elasticity
