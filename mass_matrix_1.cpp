#include "mass_matrix.hpp"

namespace Elasticity
{

template Eigen::Matrix<double, basis_size<C0Triangle<1>>, basis_size<C0Triangle<1>>>
element_mass_matrix(const C0Triangle<1> &);

} // namespace Elasticity
