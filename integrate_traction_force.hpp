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

#include "elasticity.hpp"
#include "mesh_tools.hpp"

namespace Elasticity
{

template <class Mesh>
Eigen::Matrix2Xd
integrate_traction_force(const Mesh &mesh, size_t bindex, const Eigen::Matrix2Xd &fc)
{
    const auto curve_mesh = get_boundary_mesh(mesh, bindex);
    constexpr int order = decltype(curve_mesh)::order();

    assert(static_cast<size_t>(fc.cols()) == curve_mesh.num_nodes());

    Eigen::Matrix2Xd F = Eigen::Matrix2Xd::Zero(2, fc.cols());

    size_t first_node = 0;
    constexpr auto mass_form = [](auto u, auto v) { return u * v; };

    for (const auto &el : curve_mesh.elements())
    {
        const auto mass_matrix = el.form_matrix(mass_form);
        Eigen::Matrix<double, order + 1, order + 1> M =
            Eigen::Map<const Eigen::Matrix<double, order + 1, order + 1>>(&mass_matrix(0, 0));

        F.block<2, order + 1>(0, first_node) += fc.block<2, order + 1>(0, first_node) * M;

        first_node += order;
    }

    return F;
}

} // namespace Elasticity

