#include "stiffness_matrix.hpp"

namespace Elasticity
{

template el_stiffness_type<C0Triangle<2>>
element_stiffness_matrix(const C0Triangle<2> &, double, double);

template el_stiffness_type<C0Triangle<2>> element_stiffness_matrix(
    const C0Triangle<2> &, const coeffs_type<C0Triangle<2>> &,
    const coeffs_type<C0Triangle<2>> &);

} // namespace Elasticity
