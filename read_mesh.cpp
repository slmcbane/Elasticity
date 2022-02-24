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

#include "read_mesh.hpp"

#include <cstdio>

#include "mesh_tools.hpp"

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
