#ifndef MESH_TRAITS_HPP
#define MESH_TRAITS_HPP

#include <type_traits>
#include <variant>

#include <Eigen/Core>

#include "Galerkin.hpp"
#include "TetMesh.hpp"

#ifndef ELASTICITY_MAX_ELEMENT_ORDER
#define ELASTICITY_MAX_ELEMENT_ORDER 2
#endif

namespace msh
{

template <class T>
struct IsTetMesh : public std::false_type
{
};

template <class CoordT, size_t... Others>
struct IsTetMesh<TetMesh<CoordT, Others...>> : public std::true_type
{
};

template <size_t NPF, size_t IN>
constexpr int which_poly_order() noexcept
{
    if constexpr (NPF == 0 && IN == 0)
    {
        return 1;
    }
    else if constexpr (NPF == 1 && IN == 0)
    {
        return 2;
    }
    else if constexpr (NPF == 2 && IN == 1)
    {
        return 3;
    }
    else if constexpr (NPF == 3 && IN == 3)
    {
        return 4;
    }
    else
    {
        return -1;
    }
}

template <class Mesh>
struct ElementOrder : public std::integral_constant<int, -1>
{
};

template <class CoordT, size_t MEA, size_t MNA, size_t NPF, size_t IN>
struct ElementOrder<TetMesh<CoordT, MEA, MNA, NPF, IN>>
    : public std::integral_constant<int, which_poly_order<NPF, IN>()>
{
};

template <class Mesh>
constexpr int element_order = ElementOrder<Mesh>::value;

} // namespace msh

namespace Elasticity
{

constexpr size_t max_node_adjacencies = 96;

constexpr size_t max_element_adjacencies = 16;

using Mesh1 = msh::TetMesh<double, max_element_adjacencies, max_node_adjacencies, 0, 0>;
using Mesh2 = msh::TetMesh<double, max_element_adjacencies, max_node_adjacencies, 1, 0>;
using Mesh3 = msh::TetMesh<double, max_element_adjacencies, max_node_adjacencies, 2, 1>;
using Mesh4 = msh::TetMesh<double, max_element_adjacencies, max_node_adjacencies, 3, 3>;

#if ELASTICITY_MAX_ELEMENT_ORDER > 3
using MeshVariant = std::variant<Mesh1, Mesh2, Mesh3, Mesh4>;
#elif ELASTICITY_MAX_ELEMENT_ORDER > 2
using MeshVariant = std::variant<Mesh1, Mesh2, Mesh3>;
#else
using MeshVariant = std::variant<Mesh1, Mesh2>;
#endif

template <int order>
using MeshType = std::decay_t<decltype(std::get<order - 1>(std::declval<MeshVariant>()))>;

template <int Order>
struct CurveMesh
{
    CurveMesh(Eigen::Matrix2Xd &&coords) : m_coords(std::move(coords))
    {
        std::array<size_t, Order + 1> tmp;
        size_t d = 0;
        double curr_point = 0;
        for (int i = 0; i < m_coords.cols() - 1; ++i)
        {
            Eigen::Vector2d diff = m_coords.col(i + 1) - m_coords.col(i);
            double dist = std::sqrt(diff.dot(diff));
            m_elements.push_back(Galerkin::Elements::IntervalElement<Order, double>(
                curr_point, curr_point + dist));
            for (size_t j = 0; j <= Order; ++j)
            {
                tmp[j] = d + j;
            }
            m_dofs.push_back(tmp);
            d += Order;
        }
    }

    const auto &elements() const noexcept { return m_elements; }

    const auto &dofs() const noexcept { return m_dofs; }

    size_t num_nodes() const noexcept { return m_elements.size() * Order + 1; }

    Eigen::Matrix2d get_element_endpoints(int eli) const
    {
        Eigen::Matrix2d endpoints;
        endpoints.col(0) = m_coords.col(eli);
        endpoints.col(1) = m_coords.col(eli + 1);
        return endpoints;
    }

    constexpr static int order() noexcept { return Order; }

  private:
    Eigen::Matrix2Xd m_coords;
    std::vector<Galerkin::Elements::IntervalElement<Order, double>> m_elements;
    std::vector<std::array<size_t, Order + 1>> m_dofs;
};

#if ELASTICITY_MAX_ELEMENT_ORDER > 3
using CurveMeshVariant = std::variant<CurveMesh<1>, CurveMesh<2>, CurveMesh<3>, CurveMesh<4>>;
#elif ELASTICITY_MAX_ELEMENT_ORDER > 2
using CurveMeshVariant = std::variant<CurveMesh<1>, CurveMesh<2>, CurveMesh<3>>;
#else
using CurveMeshVariant = std::variant<CurveMesh<1>, CurveMesh<2>>;
#endif

} // namespace Elasticity

#endif // MESH_TRAITS_HPP
