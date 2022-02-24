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

#ifndef CONVERT_TO_EIGEN_HPP
#define CONVERT_TO_EIGEN_HPP

#include "SymSparse.hpp"
#include <Eigen/Core>
#include <Eigen/SparseCore>

#include <tuple>

/*
 * Convert the SymSparse matrix to an Eigen sparse matrix.
 * Only the upper triangular part of the matrix is stored.
 */
template <typename T, size_t MaxPerRow, bool Check>
auto convert_to_eigen(const SymSparse::SymmetricSparseMatrix<T, MaxPerRow, Check> &A)
{
    using triplet = Eigen::Triplet<T>;
    std::vector<triplet> triplets;
    triplets.reserve(A.num_rows() * (MaxPerRow / 2));

    for (size_t m = 0; m < A.num_rows(); ++m)
    {
        for (const auto &entry : A.row(m, std::false_type{}))
        {
            auto n = get<0>(entry);
            triplets.emplace_back(m, n, std::get<1>(entry));
            if (n != m)
            {
                triplets.emplace_back(n, m, std::get<1>(entry));
            }
        }
    }
    Eigen::SparseMatrix<T> B(A.num_rows(), A.num_rows());
    B.setFromTriplets(triplets.begin(), triplets.end());
    return B;
}

#endif /* CONVERT_TO_EIGEN_HPP */

