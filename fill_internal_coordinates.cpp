/*
 * Copyright 2022 Sean McBane
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

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
