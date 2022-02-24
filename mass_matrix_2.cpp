#include "mass_matrix.hpp"

namespace Elasticity
{

template Eigen::Matrix<double, basis_size<C0Triangle<2>>, basis_size<C0Triangle<2>>>
element_mass_matrix(const C0Triangle<2> &);

} // namespace Elasticity
