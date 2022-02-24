/*
 * Copyright 2020 The University of Texas at Austin.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to
 * deal in the Software without restriction, including without limitation the
 * rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
 * sell copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
 * THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE.
 */

#ifndef MESH_DIAGNOSTICS_HPP
#define MESH_DIAGNOSTICS_HPP

#include <fmt/core.h>
#include <fmt/format.h>

#include "../TetMesh/TetMesh.hpp"
#include "mesh_traits.hpp"

#include <cstdio>
#include <iterator>
#include <type_traits>

namespace msh
{

/*
 * show_mesh_diagnostics(mesh, detail, file)
 *
 * Print a summary of mesh statistics to specified FILE* (default: stdout)
 * Builds the string in memory with fmt::memory_buffer.
 *
 * Specify level of detail = 0, 1, or 2. (Other values will be rounded to
 * nearest). Default is 0.
 */
template <class Mesh>
void show_mesh_diagnostics(const Mesh &mesh, int detail = 0, std::FILE *file = stdout);

template <class T>
struct num_face_nodes
{
};

template <class T>
struct num_intern_nodes
{
};

template <class CoordT, size_t MEA, size_t MNA, size_t NodesPerFace, size_t IN>
struct num_face_nodes<TetMesh<CoordT, MEA, MNA, NodesPerFace, IN>>
    : public std::integral_constant<size_t, NodesPerFace>
{
};

template <class CoordT, size_t MEA, size_t MNA, size_t NPF, size_t InternalNodes>
struct num_intern_nodes<TetMesh<CoordT, MEA, MNA, NPF, InternalNodes>>
    : public std::integral_constant<size_t, InternalNodes>
{
};

template <class T>
inline constexpr auto num_face_nodes_v = num_face_nodes<T>::value;

template <class T>
inline constexpr auto num_intern_nodes_v = num_intern_nodes<T>::value;

template <class Mesh>
auto mesh_summary_only(const Mesh &mesh)
{
    static_assert(IsTetMesh<Mesh>::value, "Passed Mesh type was not a msh::TetMesh");
    std::string buf;
    auto it = std::back_inserter(buf);
    fmt::format_to(
        it, "Mesh summary:\n"
            "=============\n");
    fmt::format_to(it, "  * {} nodes\n", mesh.num_nodes());
    fmt::format_to(
        it, "  * {} elements, with {} nodes per face and {} internal nodes\n",
        mesh.num_elements(), num_face_nodes_v<Mesh>, num_intern_nodes_v<Mesh>);
    fmt::format_to(
        it, "  * {} defined boundary segments; tags listed below:\n", mesh.num_boundaries());

    for (size_t i = 0; i < mesh.num_boundaries(); ++i)
    {
        fmt::format_to(
            it, "    - Boundary {} has {} faces; tags: ", i, mesh.boundary(i).faces.size());
        for (size_t j = 0; j < mesh.boundary_tags(i).size(); ++j)
        {
            if (j == mesh.boundary_tags(i).size() - 1)
            {
                fmt::format_to(it, "{}", mesh.boundary_tags(i)[j]);
            }
            else
            {
                fmt::format_to(it, "{}, ", mesh.boundary_tags(i)[j]);
            }
        }
        fmt::format_to(it, "\n");
    }
    return buf;
}

template <class Mesh>
auto mesh_medium_detail(const Mesh &mesh)
{
    static_assert(IsTetMesh<Mesh>::value, "Passed Mesh type was not a msh::TetMesh");
    auto buf = mesh_summary_only(mesh);
    auto it = std::back_inserter(buf);

    fmt::format_to(
        it, FMT_STRING("\nAverage bandwidth of a matrix from this mesh: {:f}\n"),
        mesh.average_bandwidth());

    fmt::format_to(it, FMT_STRING("Memory used for mesh (B): {:d}\n"), mesh.storage_used());

    return buf;
}

template <class Mesh>
auto mesh_high_detail(const Mesh &mesh)
{
    static_assert(IsTetMesh<Mesh>::value, "Passed Mesh type was not a msh::TetMesh");
    auto buf = mesh_medium_detail(mesh);
    auto it = std::back_inserter(buf);

    fmt::format_to(it, "\nNodes:\n");
    for (size_t i = 0; i < mesh.num_nodes(); ++i)
    {
        fmt::format_to(
            it, FMT_STRING("  {:d}: ({}, {})\n"), i, mesh.coord(i)[0], mesh.coord(i)[1]);
    }

    fmt::format_to(
        it, "\nElements (vertices given first, "
            "then face nodes in order, then internal):\n\n");
    for (size_t i = 0; i < mesh.num_elements(); ++i)
    {
        const auto &el = mesh.element(i);
        fmt::format_to(it, FMT_STRING("  {:d}: ("), i);
        fmt::format_to(
            it, "{}", fmt::join(el.control_nodes.begin(), el.control_nodes.end(), ", "));
        if (num_face_nodes_v < Mesh >> 0)
        {
            fmt::format_to(it, ", ");
            fmt::format_to(
                it, "{}", fmt::join(el.face_nodes[0].begin(), el.face_nodes[0].end(), ", "));
            fmt::format_to(it, ", ");
            fmt::format_to(
                it, "{}", fmt::join(el.face_nodes[1].begin(), el.face_nodes[1].end(), ", "));
            fmt::format_to(it, ", ");
            fmt::format_to(
                it, "{}", fmt::join(el.face_nodes[2].begin(), el.face_nodes[2].end(), ", "));
        }
        if (num_intern_nodes_v < Mesh >> 0)
        {
            fmt::format_to(it, ", ");
            fmt::format_to(
                it, "{}", fmt::join(el.internal_nodes.begin(), el.internal_nodes.end(), ", "));
        }
        fmt::format_to(it, ")\n");
    }

    fmt::format_to(it, "\nBoundaries (only nodes in order shown here):\n");
    for (size_t i = 0; i < mesh.num_boundaries(); ++i)
    {
        const auto &boundary = mesh.boundary(i);
        fmt::format_to(
            it, FMT_STRING("  {:d}: ({})\n"), i,
            fmt::join(boundary.nodes.begin(), boundary.nodes.end(), ", "));
    }
    return buf;
}

template <class Mesh>
void show_mesh_diagnostics(const Mesh &mesh, int detail, std::FILE *file)
{
    static_assert(IsTetMesh<Mesh>::value, "Passed Mesh type was not a msh::TetMesh");
    auto buf = detail < 1   ? mesh_summary_only(mesh)
               : detail < 2 ? mesh_medium_detail(mesh)
                            : mesh_high_detail(mesh);
    fmt::print(file, buf + "\n");
}

} // namespace msh

#endif /* MESH_DIAGNOSTICS_HPP */
