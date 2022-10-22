/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2022 by Oleg Rogozin
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE.md at
 * the top level directory of deal.II.
 *
 * ---------------------------------------------------------------------
 */

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/timer.h>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/affine_constraints.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_values.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>

#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_q.h>

#include <fstream>
#include <iostream>

#include "Parameters.h"

namespace BimetallicStrip
{
  using namespace dealii;


  template <int dim>
  class Solver
  {
  public:
    Solver(Parameters<dim> &parameters);
    void run();

  private:
    void setup_system();
    void assemble_system();
    void solve();
    void refine_mesh();
    void output_results(const unsigned int cycle) const;

    double compute_alpha(const Point<dim> &point) const;
    void   compute_alpha(const std::vector<Point<dim>> &points,
                         std::vector<double> &          values) const;

    Parameters<dim> &   params;
    const Problem<dim> &problem;

    Triangulation<dim> triangulation;
    DoFHandler<dim>    dof_handler;
    FESystem<dim>      fe;

    AffineConstraints<double> constraints;

    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> system_matrix;

    Vector<double> solution;
    Vector<double> system_rhs;
  };



  template <int dim>
  Solver<dim>::Solver(Parameters<dim> &parameters)
    : params(parameters)
    , problem(params.problem)
    , dof_handler(triangulation)
    , fe(FE_Q<dim>(params.fe.poly_degree), dim)
  {}



  template <int dim>
  void Solver<dim>::setup_system()
  {
    Timer timer;
    if (params.output.verbosity > 0)
      std::cout << " -- Setting up system... " << std::flush;

    dof_handler.distribute_dofs(fe);
    solution.reinit(dof_handler.n_dofs());
    system_rhs.reinit(dof_handler.n_dofs());

    constraints.clear();
    DoFTools::make_hanging_node_constraints(dof_handler, constraints);
    /*
    VectorTools::interpolate_boundary_values(dof_handler,
                                             0,
                                             Functions::ZeroFunction<dim>(dim),
                                             constraints);
                                             */
    constraints.close();

    DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern(dof_handler,
                                    dsp,
                                    constraints,
                                    /*keep_constrained_dofs = */ false);
    sparsity_pattern.copy_from(dsp);

    system_matrix.reinit(sparsity_pattern);

    if (params.output.verbosity > 0)
      std::cout << "done (" << timer.cpu_time() << "s)" << std::endl;
  }



  template <int dim>
  void Solver<dim>::assemble_system()
  {
    Timer timer;
    if (params.output.verbosity > 0)
      std::cout << " -- Assembling the system... " << std::flush;

    QGauss<dim> quadrature_formula(fe.degree + 1);

    FEValues<dim> fe_values(fe,
                            quadrature_formula,
                            update_gradients | update_quadrature_points |
                              update_JxW_values);

    const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
    const unsigned int n_q_points    = quadrature_formula.size();

    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
    Vector<double>     cell_rhs(dofs_per_cell);

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    std::vector<double> alpha_values(n_q_points);

    for (const auto &cell : dof_handler.active_cell_iterators())
      {
        cell_matrix = 0;
        cell_rhs    = 0;

        fe_values.reinit(cell);
        compute_alpha(fe_values.get_quadrature_points(), alpha_values);

        for (const unsigned int i : fe_values.dof_indices())
          {
            const unsigned int component_i =
              fe.system_to_component_index(i).first;

            for (const unsigned int j : fe_values.dof_indices())
              {
                const unsigned int component_j =
                  fe.system_to_component_index(j).first;

                for (const unsigned int q_point :
                     fe_values.quadrature_point_indices())
                  {
                    cell_matrix(i, j) +=
                      ((fe_values.shape_grad(i, q_point)[component_i] *
                        fe_values.shape_grad(j, q_point)[component_j] *
                        problem.lambda) +
                       (fe_values.shape_grad(i, q_point)[component_j] *
                        fe_values.shape_grad(j, q_point)[component_i] *
                        problem.mu) +
                       ((component_i == component_j) ?
                          (fe_values.shape_grad(i, q_point) *
                           fe_values.shape_grad(j, q_point) * problem.mu) :
                          0)) *
                      fe_values.JxW(q_point);
                  }
              }
          }

        for (const unsigned int i : fe_values.dof_indices())
          {
            const unsigned int component_i =
              fe.system_to_component_index(i).first;
            for (const unsigned int q_point :
                 fe_values.quadrature_point_indices())
              {
                cell_rhs(i) += fe_values.shape_grad(i, q_point)[component_i] *
                               problem.threeK * alpha_values[q_point] *
                               fe_values.JxW(q_point);
              }
          }

        cell->get_dof_indices(local_dof_indices);
        constraints.distribute_local_to_global(
          cell_matrix, cell_rhs, local_dof_indices, system_matrix, system_rhs);
      }

    if (params.output.verbosity > 0)
      std::cout << "done (" << timer.cpu_time() << "s)" << std::endl;
  }



  template <int dim>
  void Solver<dim>::solve()
  {
    Timer timer;
    if (params.output.verbosity > 0)
      std::cout << " -- Solving a linear system... " << std::flush;

    SolverControl solver_control(params.ls.max_iter,
                                 params.ls.tol * system_rhs.l2_norm());

    SolverCG<Vector<double>> cg(solver_control);

    PreconditionSSOR<SparseMatrix<double>> preconditioner;
    preconditioner.initialize(system_matrix, params.ls.preconditioner_relax);

    cg.solve(system_matrix, solution, system_rhs, preconditioner);

    if (params.output.verbosity > 0)
      std::cout << "done (" << timer.cpu_time() << "s)" << std::endl;

    std::cout << "     " << solver_control.last_step()
              << " CG iterations: initial residual = "
              << solver_control.initial_value() / system_rhs.l2_norm()
              << " final residual = "
              << solver_control.last_value() / system_rhs.l2_norm()
              << std::endl;

    constraints.distribute(solution);
  }



  template <int dim>
  void Solver<dim>::refine_mesh()
  {
    Timer timer;
    if (params.output.verbosity > 0)
      std::cout << " -- Refining the mesh... " << std::flush;

    Vector<float> estimated_error_per_cell(triangulation.n_active_cells());

    KellyErrorEstimator<dim>::estimate(dof_handler,
                                       QGauss<dim - 1>(fe.degree + 1),
                                       {},
                                       solution,
                                       estimated_error_per_cell);

    GridRefinement::refine_and_coarsen_fixed_fraction(triangulation,
                                                      estimated_error_per_cell,
                                                      params.mr.upper,
                                                      params.mr.lower,
                                                      params.mr.max_n_cells);

    triangulation.execute_coarsening_and_refinement();

    if (params.output.verbosity > 0)
      std::cout << "done (" << timer.cpu_time() << "s)" << std::endl;
  }



  template <int dim>
  void Solver<dim>::output_results(const unsigned int cycle) const
  {
    Timer timer;
    if (params.output.verbosity > 0)
      std::cout << " -- Writing results... " << std::flush;

    DataOut<dim> data_out;

    data_out.attach_dof_handler(dof_handler);

    std::vector<std::string> solution_names(dim, "displacement");

    std::vector<DataComponentInterpretation::DataComponentInterpretation>
      data_component_interpretation(
        dim, DataComponentInterpretation::component_is_part_of_vector);

    data_out.add_data_vector(solution,
                             solution_names,
                             DataOut<dim>::type_dof_data,
                             data_component_interpretation);

    data_out.build_patches();

    std::ofstream output("solution-" + std::to_string(cycle) + ".vtk");
    data_out.write_vtk(output);

    if (params.output.verbosity > 0)
      std::cout << "done (" << timer.cpu_time() << "s)" << std::endl;
  }



  template <int dim>
  void Solver<dim>::run()
  {
    for (unsigned int cycle = 0; cycle <= params.mr.j_max - params.mr.j_min;
         ++cycle)
      {
        std::cout << "Cycle " << cycle << ':' << std::endl;

        if (cycle == 0)
          {
            GridGenerator::hyper_rectangle(triangulation,
                                           Point<dim>(),
                                           problem.dimensions);
            triangulation.refine_global(params.mr.j_min);
          }
        else
          refine_mesh();

        std::cout << "   Number of active cells:       "
                  << triangulation.n_active_cells() << std::endl;

        setup_system();

        std::cout << "   Number of degrees of freedom: " << dof_handler.n_dofs()
                  << std::endl;

        assemble_system();
        solve();
        output_results(cycle);
      }
  }



  template <int dim>
  double Solver<dim>::compute_alpha(const Point<dim> &point) const
  {
    return point[problem.axis] < problem.dimensions[problem.axis] / 2 ?
             problem.alpha1 :
             problem.alpha2;
  }



  template <int dim>
  void Solver<dim>::compute_alpha(const std::vector<Point<dim>> &points,
                                  std::vector<double> &          values) const
  {
    AssertDimension(values.size(), points.size());

    for (unsigned int i = 0; i < points.size(); ++i)
      values[i] = compute_alpha(points[i]);
  }

} // namespace BimetallicStrip


int main(int argc, char *argv[])
{
  try
    {
      using namespace BimetallicStrip;

      constexpr int dim = 2;

      Parameters<dim> parameters;

      const std::string prm_file = argc > 1 ? argv[1] : "parameters.prm";
      std::cout << "Reading parameters... " << std::flush;
      ParameterAcceptor::initialize(prm_file);
      std::cout << "done" << std::endl;

      deallog.depth_console(parameters.output.verbosity);

      Solver<dim> solver(parameters);
      solver.run();
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;

      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }

  return 0;
}
