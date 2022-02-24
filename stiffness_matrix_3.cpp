#if ELASTICITY_MAX_ELEMENT_ORDER > 2

#include "stiffness_matrix.hpp"

namespace Elasticity
{

template el_stiffness_type<C0Triangle<3>>
element_stiffness_matrix(const C0Triangle<3> &, double, double);

template el_stiffness_type<C0Triangle<3>> element_stiffness_matrix(
    const C0Triangle<3> &, const coeffs_type<C0Triangle<3>> &,
    const coeffs_type<C0Triangle<3>> &);

} // namespace Elasticity

#endif
