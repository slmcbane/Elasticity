#include "stiffness_matrix.hpp"

namespace Elasticity
{

template el_stiffness_type<C0Triangle<1>>
element_stiffness_matrix(const C0Triangle<1> &, double, double);

template el_stiffness_type<C0Triangle<1>> element_stiffness_matrix(
    const C0Triangle<1> &, const coeffs_type<C0Triangle<1>> &,
    const coeffs_type<C0Triangle<1>> &);

} // namespace Elasticity
