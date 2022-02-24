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

#ifndef GET_TRIANGLES_HPP
#define GET_TRIANGLES_HPP

#include "mesh_traits.hpp"
#include <iterator>

namespace Elasticity
{

/*
 * Contains the most basic info extracted from a mesh - a bunch of coordinates
 * of nodes, and a bunch of triangles specified as a 3-tuple of their vertices,
 * as indexed in the array of coordinates. This is simpler than the TetMesh
 * representation and we lower to this representation to save mesh info for
 * postprocessing.
 *
 * Vertex indices are stored as signed 64-bit integers for interop with Julia.
 */
struct TriangleInfo
{
    std::vector<std::array<double, 2>> coords;
    std::vector<std::array<int64_t, 3>> vertices;

    void serialize(FILE *out) const
    {
        int64_t num_nodes = coords.size();
        size_t count = fwrite(&num_nodes, sizeof(int64_t), 1, out);
        if (count != 1)
        {
            fprintf(stderr, "Failed at serializing TriangleInfo!\n");
            return;
        }
        int64_t num_elements = vertices.size();
        count = fwrite(&num_elements, sizeof(int64_t), 1, out);
        if (count != 1)
        {
            fprintf(stderr, "Failed at serializing TriangleInfo!\n");
            return;
        }

        count = fwrite(coords.data(), sizeof(std::array<double, 2>), coords.size(), out);
        if (count != coords.size())
        {
            fprintf(stderr, "Failed at serializing TriangleInfo!\n");
            return;
        }

        count = fwrite(vertices.data(), sizeof(std::array<int64_t, 3>), num_elements, out);
        if (count != vertices.size())
        {
            fprintf(stderr, "Failed at serializing TriangleInfo!\n");
            return;
        }
    }
};

/*
 * Convert one of the 4 mesh types to TriangleInfo. This decomposes higher-order
 * elements into sub-triangles.
 */
template <class Mesh>
TriangleInfo get_triangle_info(const Mesh &mesh)
{
    constexpr int order = msh::element_order<Mesh>;
    TriangleInfo tinfo;
    tinfo.coords.reserve(mesh.num_nodes());

    for (size_t ni = 0; ni < mesh.num_nodes(); ++ni)
    {
        tinfo.coords.emplace_back(mesh.coord(ni));
    }

    auto insert_vert = [&tinfo](size_t i, size_t j, size_t k)
    {
        tinfo.vertices.emplace_back();
        tinfo.vertices.back()[0] = i;
        tinfo.vertices.back()[1] = j;
        tinfo.vertices.back()[2] = k;
    };

    if constexpr (order == 1)
    {
        tinfo.vertices.reserve(mesh.num_elements());
        for (size_t eli = 0; eli < mesh.num_elements(); ++eli)
        {
            const auto &el_info = mesh.element(eli);
            insert_vert(
                el_info.control_nodes[0], el_info.control_nodes[1], el_info.control_nodes[2]);
        }
    }
    else if constexpr (order == 2)
    {
        tinfo.vertices.reserve(4 * mesh.num_elements());
        for (size_t eli = 0; eli < mesh.num_elements(); ++eli)
        {
            const auto nn = mesh.element(eli).node_numbers();
            insert_vert(nn[0], nn[3], nn[5]);
            insert_vert(nn[3], nn[1], nn[4]);
            insert_vert(nn[4], nn[2], nn[5]);
            insert_vert(nn[5], nn[3], nn[4]);
        }
    }
    else if constexpr (order == 3)
    {
        tinfo.vertices.reserve(9 * mesh.num_elements());
        for (size_t eli = 0; eli < mesh.num_elements(); ++eli)
        {
            const auto nn = mesh.element(eli).node_numbers();
            insert_vert(nn[0], nn[3], nn[8]);
            insert_vert(nn[8], nn[3], nn[9]);
            insert_vert(nn[9], nn[7], nn[8]);
            insert_vert(nn[6], nn[7], nn[9]);
            insert_vert(nn[6], nn[2], nn[7]);
            insert_vert(nn[9], nn[3], nn[4]);
            insert_vert(nn[5], nn[9], nn[4]);
            insert_vert(nn[5], nn[6], nn[9]);
            insert_vert(nn[1], nn[5], nn[4]);
        }
    }
    else
    {
        tinfo.vertices.reserve(16 * mesh.num_elements());
        for (size_t eli = 0; eli < mesh.num_elements(); ++eli)
        {
            const auto nn = mesh.element(eli).node_numbers();
            insert_vert(nn[0], nn[3], nn[11]);
            insert_vert(nn[11], nn[3], nn[12]);
            insert_vert(nn[10], nn[11], nn[12]);
            insert_vert(nn[10], nn[12], nn[14]);
            insert_vert(nn[9], nn[10], nn[14]);
            insert_vert(nn[8], nn[9], nn[14]);
            insert_vert(nn[2], nn[9], nn[8]);
            insert_vert(nn[12], nn[3], nn[4]);
            insert_vert(nn[12], nn[4], nn[13]);
            insert_vert(nn[14], nn[12], nn[13]);
            insert_vert(nn[14], nn[13], nn[7]);
            insert_vert(nn[8], nn[14], nn[7]);
            insert_vert(nn[13], nn[4], nn[5]);
            insert_vert(nn[13], nn[5], nn[6]);
            insert_vert(nn[7], nn[13], nn[6]);
            insert_vert(nn[6], nn[5], nn[1]);
        }
    }
    return tinfo;
}

template <>
inline TriangleInfo get_triangle_info(const MeshVariant &mv)
{
    return std::visit([](const auto &mesh) { return get_triangle_info(mesh); }, mv);
}

} // namespace Elasticity

#endif // GET_TRIANGLES_HPP

