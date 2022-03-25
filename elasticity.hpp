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

#ifndef LINEAR_ELASTICITY_HPP
#define LINEAR_ELASTICITY_HPP

#include <Eigen/Core>
#include <cassert>
#include <exception>

#include "C0_triangles.hpp"
#include "Galerkin.hpp"
#include "SymSparse.hpp"
#include "mesh_traits.hpp"
#include "src/FunctionBase.hpp"

namespace Elasticity
{

template <class Element>
struct basis_size_struct
{
};

template <>
struct basis_size_struct<C0Triangle<1>>
{
    constexpr static size_t value = 3;
};

template <>
struct basis_size_struct<C0Triangle<2>>
{
    constexpr static size_t value = 6;
};

template <>
struct basis_size_struct<C0Triangle<3>>
{
    constexpr static size_t value = 10;
};

template <>
struct basis_size_struct<C0Triangle<4>>
{
    constexpr static size_t value = 15;
};

template <class Element>
constexpr size_t basis_size = basis_size_struct<Element>::value;

template <class Element>
using coeffs_type = Eigen::Matrix<double, basis_size<Element>, 1>;

template <class Element>
using vector_coeffs_type = Eigen::Matrix<double, 2, basis_size<Element>>;

template <class Mesh>
auto get_scalar_coeffs(const Mesh &mesh, size_t eli, const Eigen::VectorXd &uc)
{
    const auto &el_info = mesh.element(eli);
    Eigen::Matrix<double, std::decay_t<decltype(el_info)>::num_nodes(), 1> coeffs;

    auto nn = el_info.node_numbers();
    for (size_t i = 0; i < el_info.num_nodes(); ++i)
    {
        coeffs[i] = uc[nn[i]];
    }

    return coeffs;
}

template <class Mesh>
auto get_vector_coeffs(const Mesh &mesh, size_t eli, const Eigen::Matrix2Xd &uc)
{
    const auto &el_info = mesh.element(eli);
    Eigen::Matrix<double, 2, std::decay_t<decltype(el_info)>::num_nodes()> coeffs;

    auto nn = el_info.node_numbers();
    for (size_t i = 0; i < el_info.num_nodes(); ++i)
    {
        coeffs.col(i) = uc.col(nn[i]);
    }

    return coeffs;
}

/*
 * Given the coefficients of a function u in the finite element basis for el,
 * construct a function object that evaluates that function (and whose
derivatives,
 * etc. can be used).
 */
template <class Element>
auto make_function(const Element &, const coeffs_type<Element> &uc)
{
    using Galerkin::Functions::ConstantFunction;
    return Galerkin::static_sum<1, basis_size<Element>>(
        [&uc](auto I) { return ConstantFunction(uc[I()]) * Galerkin::get<I()>(Element::basis); },
        ConstantFunction(uc[0]) * Galerkin::get<0>(Element::basis));
}

/*
 * Given the coefficients of a *vector* function u in the finite element basis
 * for el, construct 2 function objects evaluating the X and Y parts of the
 * vector function.
 */
template <class Element>
auto make_vector_functions(const Element &el, const vector_coeffs_type<Element> &uc)
{
    using Galerkin::Functions::ConstantFunction;

    auto ux = Galerkin::static_sum<1, basis_size<Element>>(
        [&](auto I) { return ConstantFunction(uc(0, I())) * Galerkin::get<I()>(Element::basis); },
        ConstantFunction(uc(0, 0)) * Galerkin::get<0>(Element::basis));
    auto uy = Galerkin::static_sum<1, basis_size<Element>>(
        [&](auto I) { return ConstantFunction(uc(1, I())) * Galerkin::get<I()>(Element::basis); },
        ConstantFunction(uc(1, 0)) * Galerkin::get<0>(Element::basis));

    return std::pair(ux, uy);
}

template <class Element>
using el_stiffness_type = Eigen::Matrix<double, 2 * basis_size<Element>, 2 * basis_size<Element>>;

template <class Element>
el_stiffness_type<Element> element_stiffness_matrix(const Element &el, double lambda, double mu);

extern template el_stiffness_type<C0Triangle<1>>
element_stiffness_matrix(const C0Triangle<1> &, double, double);
extern template el_stiffness_type<C0Triangle<2>>
element_stiffness_matrix(const C0Triangle<2> &, double, double);

#if ELASTICITY_MAX_ELEMENT_ORDER > 2
extern template el_stiffness_type<C0Triangle<3>>
element_stiffness_matrix(const C0Triangle<3> &, double, double);
#endif

#if ELASTICITY_MAX_ELEMENT_ORDER > 3
extern template el_stiffness_type<C0Triangle<4>>
element_stiffness_matrix(const C0Triangle<4> &, double, double);
#endif

template <class Element>
el_stiffness_type<Element> element_stiffness_matrix(
    const Element &el, const coeffs_type<Element> &lambda, const coeffs_type<Element> &mu);

extern template el_stiffness_type<C0Triangle<1>> element_stiffness_matrix(
    const C0Triangle<1> &, const coeffs_type<C0Triangle<1>> &,
    const coeffs_type<C0Triangle<1>> &);
extern template el_stiffness_type<C0Triangle<2>> element_stiffness_matrix(
    const C0Triangle<2> &, const coeffs_type<C0Triangle<2>> &,
    const coeffs_type<C0Triangle<2>> &);

#if ELASTICITY_MAX_ELEMENT_ORDER > 2
extern template el_stiffness_type<C0Triangle<3>> element_stiffness_matrix(
    const C0Triangle<3> &, const coeffs_type<C0Triangle<3>> &,
    const coeffs_type<C0Triangle<3>> &);
#endif

#if ELASTICITY_MAX_ELEMENT_ORDER > 3
extern template el_stiffness_type<C0Triangle<4>> element_stiffness_matrix(
    const C0Triangle<4> &, const coeffs_type<C0Triangle<4>> &,
    const coeffs_type<C0Triangle<4>> &);
#endif

template <class Element>
Eigen::Matrix<double, basis_size<Element>, basis_size<Element>>
element_mass_matrix(const Element &el);

extern template Eigen::Matrix<double, basis_size<C0Triangle<1>>, basis_size<C0Triangle<1>>>
element_mass_matrix(const C0Triangle<1> &);
extern template Eigen::Matrix<double, basis_size<C0Triangle<2>>, basis_size<C0Triangle<2>>>
element_mass_matrix(const C0Triangle<2> &);

#if ELASTICITY_MAX_ELEMENT_ORDER > 2
extern template Eigen::Matrix<double, basis_size<C0Triangle<3>>, basis_size<C0Triangle<3>>>
element_mass_matrix(const C0Triangle<3> &);
#endif

#if ELASTICITY_MAX_ELEMENT_ORDER > 3
extern template Eigen::Matrix<double, basis_size<C0Triangle<4>>, basis_size<C0Triangle<4>>>
element_mass_matrix(const C0Triangle<4> &);
#endif

template <class Mesh>
auto instantiate_element(const Mesh &mesh, size_t which)
{
    constexpr int order = msh::element_order<Mesh>;
    const auto &el = mesh.element(which);
    return C0Triangle<order>(
        mesh.coord(el.control_nodes[0]), mesh.coord(el.control_nodes[1]),
        mesh.coord(el.control_nodes[2]));
}

constexpr auto lame_parameters(double E, double nu) noexcept
{
    double lambda = E * nu / ((1 + nu) * (1 - nu));
    double mu = E / (2 * (1 + nu));
    return std::make_pair(lambda, mu);
}

typedef SymSparse::SymmetricSparseMatrix<double, 2 * max_node_adjacencies> StiffnessType;

template <class Mesh>
StiffnessType assemble_stiffness(const Mesh &mesh, double E, double nu)
{
    StiffnessType K(mesh.num_nodes() * 2);
    auto [lambda, mu] = lame_parameters(E, nu);

    for (size_t i = 0; i < mesh.num_elements(); ++i)
    {
        const auto &el_info = mesh.element(i);
        const auto nn = el_info.node_numbers();
        const auto el = instantiate_element(mesh, i);
        const auto local_K = element_stiffness_matrix(el, lambda, mu);
        for (size_t j = 0; j < nn.size(); ++j)
        {
            for (size_t k = j; k < nn.size(); ++k)
            {
                K.insert_entry(nn[j] * 2, nn[k] * 2, local_K(2 * j, 2 * k));
                if (k != j)
                {
                    K.insert_entry(nn[j] * 2 + 1, nn[k] * 2, local_K(2 * j + 1, 2 * k));
                }
                K.insert_entry(nn[j] * 2, nn[k] * 2 + 1, local_K(2 * j, 2 * k + 1));
                K.insert_entry(nn[j] * 2 + 1, nn[k] * 2 + 1, local_K(2 * j + 1, 2 * k + 1));
            }
        }
    }
    return K;
}

template <class Mesh>
StiffnessType
assemble_stiffness(const Mesh &mesh, const Eigen::VectorXd &lambda, const Eigen::VectorXd &mu)
{
    StiffnessType K(mesh.num_nodes() * 2);

    for (size_t eli = 0; eli < mesh.num_elements(); ++eli)
    {
        auto lc = get_scalar_coeffs(mesh, eli, lambda);
        auto uc = get_scalar_coeffs(mesh, eli, mu);
        const auto &el_info = mesh.element(eli);
        const auto nn = el_info.node_numbers();
        const auto el = instantiate_element(mesh, eli);
        const auto local_K = element_stiffness_matrix(el, lc, uc);

        for (size_t j = 0; j < nn.size(); ++j)
        {
            for (size_t k = j; k < nn.size(); ++k)
            {
                K.insert_entry(nn[j] * 2, nn[k] * 2, local_K(2 * j, 2 * k));
                if (k != j)
                {
                    K.insert_entry(nn[j] * 2 + 1, nn[k] * 2, local_K(2 * j + 1, 2 * k));
                }
                K.insert_entry(nn[j] * 2, nn[k] * 2 + 1, local_K(2 * j, 2 * k + 1));
                K.insert_entry(nn[j] * 2 + 1, nn[k] * 2 + 1, local_K(2 * j + 1, 2 * k + 1));
            }
        }
    }

    return K;
}

extern template StiffnessType assemble_stiffness<Mesh1>(const Mesh1 &, double, double);
extern template StiffnessType
assemble_stiffness<Mesh1>(const Mesh1 &, const Eigen::VectorXd &, const Eigen::VectorXd &);
extern template StiffnessType assemble_stiffness<Mesh2>(const Mesh2 &, double, double);
extern template StiffnessType
assemble_stiffness<Mesh2>(const Mesh2 &, const Eigen::VectorXd &, const Eigen::VectorXd &);

#if ELASTICITY_MAX_ELEMENT_ORDER > 2
extern template StiffnessType assemble_stiffness<Mesh3>(const Mesh3 &, double, double);
extern template StiffnessType
assemble_stiffness<Mesh3>(const Mesh3 &, const Eigen::VectorXd &, const Eigen::VectorXd &);
#endif

#if ELASTICITY_MAX_ELEMENT_ORDER > 3
extern template StiffnessType assemble_stiffness<Mesh4>(const Mesh4 &, double, double);
extern template StiffnessType
assemble_stiffness<Mesh4>(const Mesh4 &, const Eigen::VectorXd &, const Eigen::VectorXd &);
#endif

template <>
inline StiffnessType
assemble_stiffness<MeshVariant>(const MeshVariant &mv, double lambda, double mu)
{
    return std::visit(
        [lambda, mu](const auto &mesh) { return assemble_stiffness(mesh, lambda, mu); }, mv);
}

template <>
inline StiffnessType assemble_stiffness<MeshVariant>(
    const MeshVariant &mv, const Eigen::VectorXd &lambda, const Eigen::VectorXd &mu)
{
    return std::visit(
        [&lambda, &mu](const auto &mesh) { return assemble_stiffness(mesh, lambda, mu); }, mv);
}

/*
 * Add a homogeneous Dirichlet condition on the boundary segment of the mesh
 * indexed by `which`. K is the stiffness matrix obtained from
 * `assemble_stiffness`, and `rhs` is the forcing vector.
 */
template <class Mesh, class RHS>
void impose_homogeneous_condition(
    const Mesh &mesh, StiffnessType &K, RHS &rhs, size_t which, double scale = 1.0)
{
    std::vector<size_t> adjacent;
    adjacent.reserve(2 * max_node_adjacencies);
    const auto &boundary = mesh.boundary(which);

    for (auto n : boundary.nodes)
    {
        adjacent.clear();
        // Get all of the adjacent DOFs to this one; since there are two
        // components of displacement there are 2 degrees of freedom (2*n,
        // 2*n+1) corresponding to each node.
        for (auto n2 : mesh.adjacent_nodes(n))
        {
            adjacent.push_back(2 * n2);
            adjacent.push_back(2 * n2 + 1);
        }
        K.eliminate_dof(2 * n, 0.0, scale, rhs, adjacent);
        K.eliminate_dof(2 * n + 1, 0.0, scale, rhs, adjacent);
    }
}

template <class Mesh>
void impose_homogeneous_condition(
    const Mesh &mesh, StiffnessType &K, size_t which, double scale = 1.0)
{
    std::vector<size_t> adjacent;
    adjacent.reserve(2 * max_node_adjacencies);
    const auto &boundary = mesh.boundary(which);

    for (auto n : boundary.nodes)
    {
        adjacent.clear();
        // Get all of the adjacent DOFs to this one; since there are two
        // components of displacement there are 2 degrees of freedom (2*n,
        // 2*n+1) corresponding to each node.
        for (auto n2 : mesh.adjacent_nodes(n))
        {
            adjacent.push_back(2 * n2);
            adjacent.push_back(2 * n2 + 1);
        }
        K.eliminate_dof(2 * n, 0.0, scale, adjacent);
        K.eliminate_dof(2 * n + 1, 0.0, scale, adjacent);
    }
}

template <>
void impose_homogeneous_condition(const MeshVariant &, StiffnessType &, size_t, double);

struct OutOfBoundsIndex : public std::exception
{
    const char *msg;
    OutOfBoundsIndex(const char *m) : msg(m) {}
    const char *what() const noexcept { return msg; }
};

template <class Mesh>
void impose_dirichlet_condition(
    const Mesh &mesh, StiffnessType &K, Eigen::Matrix2Xd &forcing, size_t which,
    const Eigen::Matrix2Xd &value, double scale = 1.0);

extern template void impose_dirichlet_condition(
    const Mesh1 &, StiffnessType &, Eigen::Matrix2Xd &, size_t, const Eigen::Matrix2Xd &, double);
extern template void impose_dirichlet_condition(
    const Mesh2 &, StiffnessType &, Eigen::Matrix2Xd &, size_t, const Eigen::Matrix2Xd &, double);

#if ELASTICITY_MAX_ELEMENT_ORDER > 2
extern template void impose_dirichlet_condition(
    const Mesh3 &, StiffnessType &, Eigen::Matrix2Xd &, size_t, const Eigen::Matrix2Xd &, double);
#endif

#if ELASTICITY_MAX_ELEMENT_ORDER > 3
extern template void impose_dirichlet_condition(
    const Mesh4 &, StiffnessType &, Eigen::Matrix2Xd &, size_t, const Eigen::Matrix2Xd &, double);
#endif

template <>
inline void impose_dirichlet_condition(
    const MeshVariant &mv, StiffnessType &K, Eigen::Matrix2Xd &forcing, size_t which,
    const Eigen::Matrix2Xd &value, double scale)
{
    std::visit(
        [&, scale, which](const auto &mesh)
        { impose_dirichlet_condition(mesh, K, forcing, which, value, scale); },
        mv);
}

template <class Force, class RHS>
void add_point_force(const Force &force, size_t node, RHS &F)
{
    if (2 * node + 1 > F.size())
    {
        throw OutOfBoundsIndex("Node is out of bounds for given forcing vector");
    }
    F[2 * node] += force[0];
    F[2 * node + 1] += force[1];
}

template <class T, class RHS, int N, class IndexContainer, int... Options>
void add_point_forces(
    const Eigen::Matrix<T, 2, N, Options...> &force, const IndexContainer &nodes, RHS &F)
{
    static_assert(sizeof...(Options) == 3);
    size_t col = 0;
    for (size_t n : nodes)
    {
        add_point_force(force.col(col), n, F);
        col += 1;
    }
}

/*
 * Where fc is a matrix containing the coefficients of the X and Y components of
 * the volume forcing in the first and second rows, respectively, integrate this
 * volume force to obtain its contribution to the linear form vector.
 */
template <class Mesh>
Eigen::Matrix2Xd integrate_volume_force(const Mesh &mesh, const Eigen::Matrix2Xd &fc);

extern template Eigen::Matrix2Xd integrate_volume_force(const Mesh1 &, const Eigen::Matrix2Xd &);
extern template Eigen::Matrix2Xd integrate_volume_force(const Mesh2 &, const Eigen::Matrix2Xd &);

#if ELASTICITY_MAX_ELEMENT_ORDER > 2
extern template Eigen::Matrix2Xd integrate_volume_force(const Mesh3 &, const Eigen::Matrix2Xd &);
#endif

#if ELASTICITY_MAX_ELEMENT_ORDER > 3
extern template Eigen::Matrix2Xd integrate_volume_force(const Mesh4 &, const Eigen::Matrix2Xd &);
#endif

template <>
inline Eigen::Matrix2Xd integrate_volume_force(const MeshVariant &mv, const Eigen::Matrix2Xd &fc)
{
    return std::visit([&fc](const auto &mesh) { return integrate_volume_force(mesh, fc); }, mv);
}

/*
 * Where fc is a matrix containing coefficients of the X and Y parts of a
 * traction force on boundary bound_index, integrates this traction force to
 * find its contribution to the linear form vector.
 */
template <class Mesh>
Eigen::Matrix2Xd
integrate_traction_force(const Mesh &mesh, size_t bound_index, const Eigen::Matrix2Xd &fc);

extern template Eigen::Matrix2Xd
integrate_traction_force(const Mesh1 &, size_t, const Eigen::Matrix2Xd &);
extern template Eigen::Matrix2Xd
integrate_traction_force(const Mesh2 &, size_t, const Eigen::Matrix2Xd &);

#if ELASTICITY_MAX_ELEMENT_ORDER > 2
extern template Eigen::Matrix2Xd
integrate_traction_force(const Mesh3 &, size_t, const Eigen::Matrix2Xd &);
#endif

#if ELASTICITY_MAX_ELEMENT_ORDER > 3
extern template Eigen::Matrix2Xd
integrate_traction_force(const Mesh4 &, size_t, const Eigen::Matrix2Xd &);
#endif

template <>
inline Eigen::Matrix2Xd
integrate_traction_force(const MeshVariant &mv, size_t bound_index, const Eigen::Matrix2Xd &fc)
{
    return std::visit(
        [&fc, bound_index](const auto &mesh)
        { return integrate_traction_force(mesh, bound_index, fc); },
        mv);
}

/*
 * Where 'traction' is an integrated traction force on boundary 'bindex',
 * as returned from 'integrate_traction_force', adds this force to the
 * appropriate columns of 'forcing' (which, for example, might have come from
 * 'integrate_volume_force'.)
 */
void add_traction_force(
    const MeshVariant &mv, Eigen::Matrix2Xd &forcing, size_t bindex,
    const Eigen::Matrix2Xd &traction);

/*
 * Where u_true is a function taking X and Y coordinates and u is a matrix with
 * coefficients of X and Y parts of a FEM solution in its columns, compute the
 * L^2 norm of u - u_true and of u_true. Used for manufactured solution
 * verification.
 */
std::pair<double, double> integrate_error(
    const MeshVariant &mesh, const Eigen::Matrix2Xd &u,
    Eigen::Vector2d (*u_true)(std::array<double, 2>));

} // namespace Elasticity

#endif // LINEAR_ELASTICITY_HPP
