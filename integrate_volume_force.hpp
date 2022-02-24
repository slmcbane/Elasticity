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

namespace Elasticity
{

template <class Mesh>
Eigen::Matrix2Xd integrate_volume_force(const Mesh &mesh, const Eigen::Matrix2Xd &fc)
{
    assert(static_cast<size_t>(fc.cols()) == mesh.num_nodes());
    Eigen::Matrix2Xd F = Eigen::Matrix2Xd::Zero(2, mesh.num_nodes());

    for (size_t eli = 0; eli < mesh.num_elements(); ++eli)
    {
        const auto el = instantiate_element(mesh, eli);
        const auto nn = mesh.element(eli).node_numbers();
        auto elfc = get_vector_coeffs(mesh, eli, fc);
        auto M = element_mass_matrix(el);
        decltype(elfc) integrated_force = elfc * M;

        assert(static_cast<size_t>(integrated_force.cols()) == nn.size());
        for (size_t i = 0; i < nn.size(); ++i)
        {
            F.col(nn[i]) += integrated_force.col(i);
        }
    }

    return F;
}

} // namespace Elasticity
