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

namespace detail
{
template <class Mesh>
void add_traction_force(
    const Mesh &mesh, Eigen::Matrix2Xd &nodal_forcing, size_t bindex,
    const Eigen::Matrix2Xd &traction)
{
    const auto &brep = mesh.boundary(bindex);
    size_t bnode = 0;
    for (size_t node : brep.nodes)
    {
        nodal_forcing.col(node) += traction.col(bnode++);
    }
}

} // namespace detail

void add_traction_force(
    const MeshVariant &mv, Eigen::Matrix2Xd &nodal_forcing, size_t bindex,
    const Eigen::Matrix2Xd &traction)
{
    std::visit(
        [&, bindex](const auto &mesh)
        { detail::add_traction_force(mesh, nodal_forcing, bindex, traction); },
        mv);
}

} // namespace Elasticity

