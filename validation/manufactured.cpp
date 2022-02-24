/*
 * Copyright 2022 Sean McBane under the terms of the MIT license.
 * This file contains validation via method of manufactured solution for the
 * linear elasticity code. The analytical solution for displacement is given
 * by
 *
 * u_x = (-3.0/2.0*pow(x, 3) + (23.0/4.0)*pow(x, 2) - 5*x)*sin(M_PI*y) +
 *  ((3.0/2.0)*pow(x, 3) - 13.0/4.0*pow(x, 2) + 1)*cos(M_PI*y);
 *
 * u_y = -2*x*y + x + 2*y - 1;
 *
 * Young's modulus: E = 1 + sin(M_PI * x) * cos(M_PI * y) / 4;
 * Poisson's ratio: nu = 3 / 10;
 *
 * The functions below are all derived symbolically from these solutions.
 */

#include <Eigen/Core>
#include <Eigen/SparseCholesky>
#include <cassert>
#include <chrono>
#include <cmath>

#include "convert_to_eigen.hpp"
#include "elasticity.hpp"
#include "mesh_tools.hpp"
#include "read_mesh.hpp"

using namespace Elasticity;
using namespace std::chrono;

double manufactured_ux(double x, double y)
{
    return (-3.0 / 2.0 * pow(x, 3) + (23.0 / 4.0) * pow(x, 2) - 5 * x) * sin(M_PI * y) +
           ((3.0 / 2.0) * pow(x, 3) - 13.0 / 4.0 * pow(x, 2) + 1) * cos(M_PI * y);
}

double manufactured_uy(double x, double y) { return -2 * x * y + x + 2 * y - 1; }

Eigen::Vector2d u_true(std::array<double, 2> coord)
{
    const auto [x, y] = coord;
    return Eigen::Vector2d(manufactured_ux(x, y), manufactured_uy(x, y));
}

Eigen::Vector2d u_right(std::array<double, 2> coord) { return u_true(coord); }

Eigen::Vector2d u_left(std::array<double, 2> coord) { return u_true(coord); }

Eigen::Vector2d u_upper(std::array<double, 2> coord) { return u_true(coord); }

Eigen::Vector2d u_lower(std::array<double, 2> coord) { return u_true(coord); }

/*
 * Definitions of the Lame parameters for the manufactured solution.
 */
double manufactured_lambda(std::array<double, 2> coord)
{
    auto [x, y] = coord;
    return (15.0 / 182.0) * sin(M_PI * x) * cos(M_PI * y) + 30.0 / 91.0;
}

double manufactured_mu(std::array<double, 2> coord)
{
    auto [x, y] = coord;
    return (5.0 / 52.0) * sin(M_PI * x) * cos(M_PI * y) + 5.0 / 13.0;
}

/*
 * Manufactured solution traction boundary condition for the lower boundary.
 */
double traction_lower_x(double x, [[maybe_unused]] double y) noexcept
{
    return (5.0 / 208.0) * (M_PI * x * (6 * pow(x, 2) - 23 * x + 20) - 4) * (sin(M_PI * x) + 4);
}

double traction_lower_y(double x, [[maybe_unused]] double y) noexcept
{
    return (5.0 / 364.0) * (sin(M_PI * x) + 4) * (-27 * pow(x, 2) + 79 * x - 40);
}

Eigen::Vector2d traction_lower(std::array<double, 2> coord) noexcept
{
    auto [x, y] = coord;
    return Eigen::Vector2d(traction_lower_x(x, y), traction_lower_y(x, y));
}

/*
 * Manufactured solution traction boundary condition for the upper boundary.
 */
double traction_upper_x(double x, [[maybe_unused]] double y) noexcept
{
    return (5.0 / 208.0) * (-M_PI * x * (6 * pow(x, 2) - 23 * x + 20) + 4) * (sin(M_PI * x) - 4);
}

double traction_upper_y(double x, [[maybe_unused]] double y) noexcept
{
    return (5.0 / 364.0) * (sin(M_PI * x) - 4) * (27 * pow(x, 2) + x - 40);
}

Eigen::Vector2d traction_upper(std::array<double, 2> coord) noexcept
{
    auto [x, y] = coord;
    return Eigen::Vector2d(traction_upper_x(x, y), traction_upper_y(x, y));
}

/*
 * Manufactured solution forcing term on \Omega
 */
double forcing_x(double x, double y) noexcept
{
    return (5.0 / 52.0) * ((13 - 18 * x) * cos(M_PI * y) + (18 * x - 23) * sin(M_PI * y)) *
               (sin(M_PI * x) * cos(M_PI * y) + 4) +
           (15.0 / 364.0) * (sin(M_PI * x) * cos(M_PI * y) + 4) *
               ((13 - 18 * x) * cos(M_PI * y) + (18 * x - 23) * sin(M_PI * y) + 4) +
           (5.0 / 208.0) * (sin(M_PI * x) * cos(M_PI * y) + 4) *
               (-pow(M_PI, 2) * x * (6 * pow(x, 2) - 23 * x + 20) * sin(M_PI * y) +
                pow(M_PI, 2) * (6 * pow(x, 3) - 13 * pow(x, 2) + 4) * cos(M_PI * y) + 8) -
           5.0 / 52.0 * M_PI *
               (x * (9 * x - 13) * cos(M_PI * y) +
                (-9 * pow(x, 2) + 23 * x - 10) * sin(M_PI * y)) *
               cos(M_PI * x) * cos(M_PI * y) +
           (15.0 / 364.0) * M_PI *
               (-x * (9 * x - 13) * cos(M_PI * y) + 4 * x +
                (9 * pow(x, 2) - 23 * x + 10) * sin(M_PI * y) - 4) *
               cos(M_PI * x) * cos(M_PI * y) -
           5.0 / 208.0 * M_PI *
               (M_PI * x * (6 * pow(x, 2) - 23 * x + 20) * cos(M_PI * y) + 8 * y +
                M_PI * (6 * pow(x, 3) - 13 * pow(x, 2) + 4) * sin(M_PI * y) - 4) *
               sin(M_PI * x) * sin(M_PI * y);
}

double forcing_y(double x, double y) noexcept
{
    return (5.0 / 1456.0) * M_PI *
           ((112 - 112 * x) * sin(M_PI * x) * sin(M_PI * y) +
            26 * (sin(M_PI * x) * cos(M_PI * y) + 4) *
                (x * (9 * x - 13) * sin(M_PI * y) +
                 (9 * pow(x, 2) - 23 * x + 10) * cos(M_PI * y)) -
            12 *
                (-x * (9 * x - 13) * cos(M_PI * y) + 4 * x +
                 (9 * pow(x, 2) - 23 * x + 10) * sin(M_PI * y) - 4) *
                sin(M_PI * x) * sin(M_PI * y) +
            7 *
                (M_PI * x * (6 * pow(x, 2) - 23 * x + 20) * cos(M_PI * y) + 8 * y +
                 M_PI * (6 * pow(x, 3) - 13 * pow(x, 2) + 4) * sin(M_PI * y) - 4) *
                cos(M_PI * x) * cos(M_PI * y));
}

Eigen::Vector2d forcing(std::array<double, 2> coord) noexcept
{
    auto [x, y] = coord;
    return Eigen::Vector2d(forcing_x(x, y), forcing_y(x, y));
}

auto compute_parameter_functions(const MeshVariant &mv)
{
    Eigen::VectorXd lambda = evaluate_on_mesh(mv, manufactured_lambda);
    Eigen::VectorXd mu = evaluate_on_mesh(mv, manufactured_mu);
    return std::pair(lambda, mu);
}

auto get_boundary_conditions(const MeshVariant &mv, size_t right_bound, size_t left_bound)
{
    Eigen::Matrix2Xd uright = evaluate_on_boundary(mv, right_bound, u_right);
    Eigen::Matrix2Xd uleft = evaluate_on_boundary(mv, left_bound, u_left);

    return std::pair(uright, uleft);
}

namespace
{

Eigen::Matrix2Xd solve(const MeshVariant &mv)
{
    size_t right_bound = find_boundary_with_tag(mv, "RIGHT").value();
    size_t left_bound = find_boundary_with_tag(mv, "LEFT").value();
    size_t upper_bound = find_boundary_with_tag(mv, "UPPER").value();
    size_t lower_bound = find_boundary_with_tag(mv, "LOWER").value();

    auto tp_start = steady_clock::now();
    Eigen::Matrix2Xd force_coefficients = evaluate_on_mesh(mv, forcing);
    Eigen::Matrix2Xd nodal_forcing = Elasticity::integrate_volume_force(mv, force_coefficients);

    Eigen::Matrix2Xd traction_upper_coefficients =
        evaluate_on_boundary(mv, upper_bound, traction_upper);
    Eigen::Matrix2Xd traction_lower_coefficients =
        evaluate_on_boundary(mv, lower_bound, traction_lower);

    Eigen::Matrix2Xd traction_upper =
        integrate_traction_force(mv, upper_bound, traction_upper_coefficients);
    Eigen::Matrix2Xd traction_lower =
        integrate_traction_force(mv, lower_bound, traction_lower_coefficients);

    add_traction_force(mv, nodal_forcing, upper_bound, traction_upper);
    add_traction_force(mv, nodal_forcing, lower_bound, traction_lower);

    auto tp_end = steady_clock::now();

    size_t us = duration_cast<microseconds>(tp_end - tp_start).count();
    printf("Assembling forcing took %ld us\n", us);

    Eigen::VectorXd lambda, mu;
    std::tie(lambda, mu) = compute_parameter_functions(mv);

    Eigen::Matrix2Xd uright, uleft;
    std::tie(uright, uleft) = get_boundary_conditions(mv, right_bound, left_bound);

    tp_start = steady_clock::now();
    auto K_unassembled = Elasticity::assemble_stiffness(mv, lambda, mu);

    Elasticity::impose_dirichlet_condition(mv, K_unassembled, nodal_forcing, right_bound, uright);
    Elasticity::impose_dirichlet_condition(mv, K_unassembled, nodal_forcing, left_bound, uleft);

    auto K = convert_to_eigen(K_unassembled);
    tp_end = steady_clock::now();
    us = duration_cast<microseconds>(tp_end - tp_start).count();
    printf("Assembling stiffness took %ld us\n", us);

    tp_start = steady_clock::now();
    Eigen::Matrix2Xd u(2, nodal_forcing.cols());
    u.reshaped() = Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>, Eigen::Upper>(K).solve(
        nodal_forcing.reshaped());
    tp_end = steady_clock::now();
    us = duration_cast<microseconds>(tp_end - tp_start).count();
    printf("Linear system solve took %ld us\n", us);

    return u;
}

} // namespace

int main(int argc, char *argv[])
{
    if (!(argc > 1))
    {
        fprintf(stderr, "Usage: ./manufactured [mesh order] [mesh files...]\n");
        return 1;
    }

    int order = atoi(argv[1]);
    if (!(order >= 1 && order <= 4))
    {
        fprintf(stderr, "Got mesh order %d; order should be 1-4\n", order);
        return 2;
    }

    for (int fi = 2; fi < argc; ++fi)
    {
        FILE *infile = fopen(argv[fi], "r");
        if (!infile)
        {
            fprintf(stderr, "Skipping mesh file %s (couldn't be opened)\n", argv[fi]);
            continue;
        }
        fclose(infile);
        auto mesh = read_mesh(argv[fi], order);
        Eigen::Matrix2Xd u = solve(mesh);

        auto [l2_error, u_norm] = Elasticity::integrate_error(mesh, u, u_true);

        printf(
            "Error for input file %s:  %E (L^2 error)  %E (L^2 u)\n", argv[fi], l2_error, u_norm);
    }

    return 0;
}

