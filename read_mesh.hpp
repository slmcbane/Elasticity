#ifndef READ_MESH_HPP
#define READ_MESH_HPP

#include <cassert>
#include <variant>

#include "mesh_traits.hpp"

namespace Elasticity
{

MeshVariant read_mesh(const char *name, int order);

} // namespace Elasticity

#endif // READ_MESH_HPP
