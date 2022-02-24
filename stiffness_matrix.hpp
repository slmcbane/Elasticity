#include "elasticity.hpp"

namespace Elasticity
{

namespace detail
{

/*
 * The contribution to the stiffness matrix from function u and v where the Lame
 * parameters are constant.
 */
template <class Element, class U, class V>
Eigen::Matrix2d
pair_stiffness_integral(const Element &el, const U &u, const V &v, double lambda, double mu)
{
    auto uxvx = el.integrate(el.template partial<0>(u) * el.template partial<0>(v));
    auto uxvy = el.integrate(el.template partial<0>(u) * el.template partial<1>(v));
    auto uyvy = el.integrate(el.template partial<1>(u) * el.template partial<1>(v));
    auto uyvx = el.integrate(el.template partial<1>(u) * el.template partial<0>(v));

    Eigen::Matrix2d K;
    K(0, 0) = (lambda + 2 * mu) * uxvx + mu * uyvy;
    K(1, 0) = lambda * uyvx + mu * uxvy;
    K(0, 1) = lambda * uxvy + mu * uyvx;
    K(1, 1) = (lambda + 2 * mu) * uyvy + mu * uxvx;

    return K;
}

/*
 * Same as above, but lambda and mu are represented by functions in the FEM
 * function space. I don't make allowance for the case where they are arbitrary
 * functions.
 */
template <class Element, class U, class V, class L, class M>
Eigen::Matrix2d
pair_stiffness_integral(const Element &el, const U &u, const V &v, const L &lambda, const M &mu)
{
    using Galerkin::Functions::ConstantFunction;

    auto uxvx = el.template partial<0>(u) * el.template partial<0>(v);
    auto uxvy = el.template partial<0>(u) * el.template partial<1>(v);
    auto uyvx = el.template partial<1>(u) * el.template partial<0>(v);
    auto uyvy = el.template partial<1>(u) * el.template partial<1>(v);

    Eigen::Matrix2d K;
    K(0, 0) = el.integrate((lambda + ConstantFunction(2) * mu) * uxvx + mu * uyvy);
    K(1, 0) = el.integrate(lambda * uyvx + mu * uxvy);
    K(0, 1) = el.integrate(lambda * uxvy + mu * uyvx);
    K(1, 1) = el.integrate((lambda + ConstantFunction(2) * mu) * uyvy + mu * uxvx);

    return K;
}

template <auto I, auto J, class Element>
void add_contribution(
    const Element &el, double lambda, double mu,
    Eigen::Matrix<double, 2 * basis_size<Element>, 2 * basis_size<Element>> &K)
{
    constexpr auto u = get<I>(Element::basis);
    constexpr auto v = get<J>(Element::basis);
    auto local_K = pair_stiffness_integral(el, u, v, lambda, mu);
    auto row = I * 2;
    auto col = J * 2;
    K.block(row, col, 2, 2) = local_K;
    K.block(col, row, 2, 2) = local_K.transpose();
}

template <auto I, auto J, class Element, class L, class M>
void add_contribution(
    const Element &el, const L &lambda, const M &mu,
    Eigen::Matrix<double, 2 * basis_size<Element>, 2 * basis_size<Element>> &K)
{
    constexpr auto u = get<I>(Element::basis);
    constexpr auto v = get<J>(Element::basis);
    auto local_K = pair_stiffness_integral(el, u, v, lambda, mu);
    auto row = I * 2;
    auto col = J * 2;
    K.block(row, col, 2, 2) = local_K;
    K.block(col, row, 2, 2) = local_K.transpose();
}

} // namespace detail

template <class Element>
el_stiffness_type<Element> element_stiffness_matrix(const Element &el, double lambda, double mu)
{
    using namespace detail;

    Eigen::Matrix<double, 2 * basis_size<Element>, 2 * basis_size<Element>> K;
    Galerkin::static_for<0, basis_size<Element>, 1>(
        [&](auto I)
        {
            Galerkin::static_for<I(), basis_size<Element>, 1>(
                [&](auto J) { add_contribution<I(), J()>(el, lambda, mu, K); });
        });
    return K;
}

template <class Element>
el_stiffness_type<Element> element_stiffness_matrix(
    const Element &el, const coeffs_type<Element> &lambda, const coeffs_type<Element> &mu)
{
    using namespace detail;

    const auto L = make_function(el, lambda);
    const auto M = make_function(el, mu);

    Eigen::Matrix<double, 2 * basis_size<Element>, 2 * basis_size<Element>> K;
    Galerkin::static_for<0, basis_size<Element>, 1>(
        [&](auto I)
        {
            Galerkin::static_for<I(), basis_size<Element>, 1>(
                [&](auto J) { add_contribution<I(), J()>(el, L, M, K); });
        });
    return K;
}

} // namespace Elasticity
