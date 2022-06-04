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

#ifndef MESH_TOOLS_HPP
#define MESH_TOOLS_HPP

#include <functional>
#include <optional>

#include "mesh_traits.hpp"

namespace Elasticity
{

template <class Mesh>
std::optional<size_t> find_boundary_with_tag(const Mesh &mesh, std::string_view tag)
{
    for (size_t i = 0; i < mesh.num_boundaries(); ++i)
    {
        const auto &tags = mesh.boundary_tags(i);
        if (std::find(tags.begin(), tags.end(), tag) != tags.end())
        {
            return i;
        }
    }
    return std::nullopt;
}

template <>
std::optional<size_t> find_boundary_with_tag(const MeshVariant &mesh, std::string_view tag);

/*
 * Get a CurveMesh of the boundary with index "bindex". Note that the curve mesh
 * is represented as flattened, with one X coordinate running from 0 to the
 * length of the boundary. Coordinates of the endpoints of each line element in
 * 2 dimensional space can be retrieved using get_element_endpoints(eli) and
 * used to do calculations in 2 dimensions.
 */
template <class Mesh>
CurveMesh<msh::element_order<Mesh>> get_boundary_mesh(const Mesh &mesh, size_t bindex)
{
    const auto &boundary = mesh.boundary(bindex);
    Eigen::Matrix2Xd coords(2, boundary.faces.size() + 1);

    auto to_eigen_vector = [](auto pt) { return Eigen::Vector2d(pt[0], pt[1]); };

    coords.col(0) = to_eigen_vector(mesh.coord(boundary.nodes[boundary.faces[0].nodes[0]]));
    int i = 1;
    for (const auto &face : boundary.faces)
    {
        size_t second_node = boundary.nodes[*(face.nodes.end() - 1)];
        coords.col(i++) = to_eigen_vector(mesh.coord(second_node));
    }

    return CurveMesh<msh::element_order<Mesh>>(std::move(coords));
}

template <class Mesh, class F>
auto evaluate_on_boundary(const Mesh &mesh, size_t which, const F &func)
{
    static_assert(
        std::is_invocable_v<F, decltype(mesh.coord(0))>,
        "func should be an object callable with a mesh coordinate");

    using result_type = decltype(func(mesh.coord(0)));
    static_assert(
        std::is_same_v<result_type, double> ||
            std::is_base_of_v<Eigen::MatrixBase<result_type>, result_type>,
        "The result type of invoking func should either be a double or an "
        "Eigen vector");

    if constexpr (std::is_same_v<result_type, double>)
    {
        const auto brep = mesh.boundary(which);
        Eigen::VectorXd f(brep.nodes.size());

        for (size_t ni = 0; ni < brep.nodes.size(); ++ni)
        {
            f[ni] = func(mesh.coord(brep.nodes[ni]));
        }
        return f;
    }
    else
    {
        constexpr int M = result_type::RowsAtCompileTime;
        static_assert(M > 0, "Size of the result must be known at compile time");

        static_assert(
            result_type::ColsAtCompileTime == 1, "Matrix-valued results are not supported");

        const auto brep = mesh.boundary(which);
        Eigen::Matrix<double, M, Eigen::Dynamic> f(M, brep.nodes.size());
        for (size_t ni = 0; ni < brep.nodes.size(); ++ni)
        {
            f.col(ni) = func(mesh.coord(brep.nodes[ni]));
        }
        return f;
    }
}

template <class F>
auto evaluate_on_boundary(const MeshVariant &mv, size_t which, const F &func)
{
    return std::visit(
        [&](const auto &mesh)
        { return evaluate_on_boundary(mesh, which, std::forward<F>(func)); },
        mv);
}

template <class Mesh, class F>
auto evaluate_on_mesh(const Mesh &mesh, const F &func)
{
    static_assert(
        std::is_invocable_v<F, decltype(mesh.coord(0))>,
        "func should be an object callable with a mesh coordinate");

    using result_type = decltype(func(mesh.coord(0)));
    static_assert(
        std::is_same_v<result_type, double> ||
            std::is_base_of_v<Eigen::MatrixBase<result_type>, result_type>,
        "The result type of invoking func should either be a double or an "
        "Eigen vector");

    if constexpr (std::is_same_v<result_type, double>)
    {
        Eigen::VectorXd result(mesh.num_nodes());
        for (size_t ni = 0; ni < mesh.num_nodes(); ++ni)
        {
            result[ni] = func(mesh.coord(ni));
        }

        return result;
    }
    else
    {
        constexpr int M = result_type::RowsAtCompileTime;
        static_assert(M > 0, "Size of the result must be known at compile time");

        static_assert(
            result_type::ColsAtCompileTime == 1, "Matrix-valued results are not supported");

        Eigen::Matrix<double, M, Eigen::Dynamic> result(M, mesh.num_nodes());

        for (size_t ni = 0; ni < mesh.num_nodes(); ++ni)
        {
            result.col(ni) = func(mesh.coord(ni));
        }
        return result;
    }
}

template <class F>
auto evaluate_on_mesh(const MeshVariant &mv, const F &func)
{
    return std::visit([&](const auto &mesh) { return evaluate_on_mesh(mesh, func); }, mv);
}

inline auto num_nodes(const MeshVariant &mv)
{
    return std::visit([](const auto &mesh) { return mesh.num_nodes(); }, mv);
}

/*
 * Higher-order meshes have nodes internal to elements to which we cannot
 * meaningfully assign coordinates in the mesh package. This function does the
 * work for third and fourth order meshes to fill the coordinates.
 */
template <class Mesh>
void fill_internal_coordinates(Mesh &mesh)
{
}

#if ELASTICITY_MAX_ELEMENT_ORDER > 2
template <>
void fill_internal_coordinates(MeshType<3> &mesh);
#endif

#if ELASTICITY_MAX_ELEMENT_ORDER > 3
template <>
void fill_internal_coordinates(MeshType<4> &mesh);
#endif

template <>
void fill_internal_coordinates(MeshVariant &mesh);

} // namespace Elasticity

#endif // MESH_TOOLS_HPP

