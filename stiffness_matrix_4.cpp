#if ELASTICITY_MAX_ELEMENT_ORDER > 3

#include "stiffness_matrix.hpp"

namespace Elasticity
{

template el_stiffness_type<C0Triangle<4>>
element_stiffness_matrix(const C0Triangle<4> &, double, double);

template el_stiffness_type<C0Triangle<4>> element_stiffness_matrix(
    const C0Triangle<4> &, const coeffs_type<C0Triangle<4>> &,
    const coeffs_type<C0Triangle<4>> &);

} // namespace Elasticity

#endif
