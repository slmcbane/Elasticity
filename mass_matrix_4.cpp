#if ELASTICITY_MAX_ELEMENT_ORDER > 3

#include "mass_matrix.hpp"

namespace Elasticity
{

template Eigen::Matrix<double, basis_size<C0Triangle<4>>, basis_size<C0Triangle<4>>>
element_mass_matrix(const C0Triangle<4> &);

} // namespace Elasticity

#endif
