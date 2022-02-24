#include "elasticity.hpp"

namespace Elasticity
{

template <class Element>
Eigen::Matrix<double, basis_size<Element>, basis_size<Element>>
element_mass_matrix(const Element &el)
{
    constexpr auto mass_form = [](const auto &u, const auto &v) { return u * v; };
    const auto form_matrix = el.form_matrix(mass_form);

    Eigen::Matrix<double, basis_size<Element>, basis_size<Element>> M;
    M = Eigen::Map<const decltype(M)>(&form_matrix(0, 0), M.rows(), M.cols());
    return M;
}

} // namespace Elasticity
