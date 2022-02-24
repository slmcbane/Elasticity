#include "elasticity.hpp"

namespace Elasticity
{

template <class Mesh>
Eigen::Matrix2Xd integrate_volume_force(const Mesh &mesh, const Eigen::Matrix2Xd &fc)
{
    assert(static_cast<size_t>(fc.cols()) == mesh.num_nodes());
    Eigen::Matrix2Xd F = Eigen::Matrix2Xd::Zero(2, mesh.num_nodes());

    for (size_t eli = 0; eli < mesh.num_elements(); ++eli)
    {
        const auto el = instantiate_element(mesh, eli);
        const auto nn = mesh.element(eli).node_numbers();
        auto elfc = get_vector_coeffs(mesh, eli, fc);
        auto M = element_mass_matrix(el);
        decltype(elfc) integrated_force = elfc * M;

        assert(static_cast<size_t>(integrated_force.cols()) == nn.size());
        for (size_t i = 0; i < nn.size(); ++i)
        {
            F.col(nn[i]) += integrated_force.col(i);
        }
    }

    return F;
}

} // namespace Elasticity
