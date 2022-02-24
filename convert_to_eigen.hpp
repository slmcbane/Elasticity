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

