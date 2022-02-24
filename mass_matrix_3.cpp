#if ELASTICITY_MAX_ELEMENT_ORDER > 2

#include "mass_matrix.hpp"

namespace Elasticity
{

template Eigen::Matrix<double, basis_size<C0Triangle<3>>, basis_size<C0Triangle<3>>>
element_mass_matrix(const C0Triangle<3> &);

} // namespace Elasticity

#endif
