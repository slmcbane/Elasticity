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
void impose_dirichlet_condition(
    const Mesh &mesh, StiffnessType &K, Eigen::Matrix2Xd &forcing, size_t which,
    const Eigen::Matrix2Xd &value, double scale)
{
    std::vector<size_t> adjacent;
    adjacent.reserve(2 * max_node_adjacencies);
    const auto &boundary = mesh.boundary(which);

    if (static_cast<size_t>(value.size()) != 2 * boundary.nodes.size())
    {
        throw OutOfBoundsIndex("Vector given for Dirichlet condition has wrong number of values");
    }

    size_t i = 0;
    auto flat_forcing = forcing.reshaped();
    for (auto n : boundary.nodes)
    {
        adjacent.clear();
        for (auto n2 : mesh.adjacent_nodes(n))
        {
            adjacent.push_back(2 * n2);
            adjacent.push_back(2 * n2 + 1);
        }
        K.eliminate_dof(2 * n, value(0, i), scale, flat_forcing, adjacent);
        K.eliminate_dof(2 * n + 1, value(1, i), scale, flat_forcing, adjacent);
        i += 1;
    }
}

} // namespace Elasticity
