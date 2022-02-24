#include "elasticity.hpp"
#include "mesh_tools.hpp"

namespace Elasticity
{

template <class Mesh>
Eigen::Matrix2Xd
integrate_traction_force(const Mesh &mesh, size_t bindex, const Eigen::Matrix2Xd &fc)
{
    const auto curve_mesh = get_boundary_mesh(mesh, bindex);
    constexpr int order = decltype(curve_mesh)::order();

    assert(static_cast<size_t>(fc.cols()) == curve_mesh.num_nodes());

    Eigen::Matrix2Xd F = Eigen::Matrix2Xd::Zero(2, fc.cols());

    size_t first_node = 0;
    constexpr auto mass_form = [](auto u, auto v) { return u * v; };

    for (const auto &el : curve_mesh.elements())
    {
        const auto mass_matrix = el.form_matrix(mass_form);
        Eigen::Matrix<double, order + 1, order + 1> M =
            Eigen::Map<const Eigen::Matrix<double, order + 1, order + 1>>(&mass_matrix(0, 0));

        F.block<2, order + 1>(0, first_node) += fc.block<2, order + 1>(0, first_node) * M;

        first_node += order;
    }

    return F;
}

} // namespace Elasticity

