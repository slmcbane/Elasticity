#include <cassert>

#include "elasticity.hpp"
#include "mesh_tools.hpp"

using namespace Galerkin;

namespace Elasticity
{

std::pair<double, double> integrate_error(
    const MeshVariant &mv, const Eigen::Matrix2Xd &u,
    Eigen::Vector2d (*u_true)(std::array<double, 2>))
{
    return std::visit(
        [&](const auto &mesh)
        {
            assert(static_cast<size_t>(u.cols()) == mesh.num_nodes());

            double l2err = 0, l2u = 0;
            const auto rule = Elasticity::instantiate_element(mesh, 0)
                                  .coordinate_map()
                                  .template quadrature_rule<8>();
            for (size_t eli = 0; eli < mesh.num_elements(); ++eli)
            {
                const auto el = Elasticity::instantiate_element(mesh, eli);
                auto ucoeffs = Elasticity::get_vector_coeffs(mesh, eli, u);
                auto [ux, uy] = Elasticity::make_vector_functions(el, ucoeffs);

                auto detJ = el.coordinate_map().detJ();

                for (size_t i = 0; i < rule.points.size(); ++i)
                {
                    Eigen::Vector2d u(ux(rule.points[i]), uy(rule.points[i]));
                    auto [x, y] = el.coordinate_map()(rule.points[i]);
                    Eigen::Vector2d utrue = u_true(std::array{x, y});
                    l2err += (u - utrue).dot(u - utrue) * rule.weights[i] * detJ(rule.points[i]);
                    l2u += utrue.dot(utrue) * rule.weights[i] * detJ(rule.points[i]);
                }
            }
            return std::pair(l2err, l2u);
        },
        mv);
}

} // namespace Elasticity

