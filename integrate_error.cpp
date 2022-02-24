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

#include <cassert>

#include "elasticity.hpp"
#include "mesh_tools.hpp"

using namespace Galerkin;

namespace Elasticity
{

std::pair<double, double> integrate_error(
    const MeshVariant &mv, const Eigen::Matrix2Xd &u,
    Eigen::Vector2d (*u_true)(std::array<double, 2>))
{
    return std::visit(
        [&](const auto &mesh)
        {
            assert(static_cast<size_t>(u.cols()) == mesh.num_nodes());

            double l2err = 0, l2u = 0;
            const auto rule = Elasticity::instantiate_element(mesh, 0)
                                  .coordinate_map()
                                  .template quadrature_rule<8>();
            for (size_t eli = 0; eli < mesh.num_elements(); ++eli)
            {
                const auto el = Elasticity::instantiate_element(mesh, eli);
                auto ucoeffs = Elasticity::get_vector_coeffs(mesh, eli, u);
                auto [ux, uy] = Elasticity::make_vector_functions(el, ucoeffs);

                auto detJ = el.coordinate_map().detJ();

                for (size_t i = 0; i < rule.points.size(); ++i)
                {
                    Eigen::Vector2d u(ux(rule.points[i]), uy(rule.points[i]));
                    auto [x, y] = el.coordinate_map()(rule.points[i]);
                    Eigen::Vector2d utrue = u_true(std::array{x, y});
                    l2err += (u - utrue).dot(u - utrue) * rule.weights[i] * detJ(rule.points[i]);
                    l2u += utrue.dot(utrue) * rule.weights[i] * detJ(rule.points[i]);
                }
            }
            return std::pair(l2err, l2u);
        },
        mv);
}

} // namespace Elasticity

