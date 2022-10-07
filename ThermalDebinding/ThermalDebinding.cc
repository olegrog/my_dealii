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

#include <deal.II/base/utilities.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/affine_constraints.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/solution_transfer.h>
#include <deal.II/numerics/matrix_tools.h>

#include <fstream>
#include <iostream>

#include "Parameters.h"

namespace ThermalDebinding
{
  using namespace dealii;


  template <int dim>
  class Solver
  {
  public:
    Solver(Parameters &parameters);
    void run();

  private:
    void   setup_system();
    void   assemble_rhs();
    void   assemble_matrix();
    double solve_ls();
    void   output_results() const;
    void   make_mesh();
    void   refine_mesh();

    // A collection of the parameters used to describe the problem setup
    Parameters &       params;
    const Problem &    problem;
    const Material &   material;
    Time &             time;
    Triangulation<dim> triangulation;
    FE_Q<dim>          fe;
    DoFHandler<dim>    dof_handler;

    AffineConstraints<double> constraints;

    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> mass_matrix;
    SparseMatrix<double> matrix;
    SparseMatrix<double> system_matrix;

    Vector<double> solution;
    Vector<double> old_rhs;
    Vector<double> rhs;
    Vector<double> system_rhs;

    double T;            // current temperature [K]
    double max_pressure; // maximum pressure in the domain [Pa]
  };



  template <int dim>
  Solver<dim>::Solver(Parameters &parameters)
    : params(parameters)
    , problem(params.problem)
    , material(params.material)
    , time(params.time)
    , fe(params.fe.poly_degree)
    , dof_handler(triangulation)
  {}



  template <int dim>
  void Solver<dim>::setup_system()
  {
    dof_handler.distribute_dofs(fe);

    std::cout << std::endl
              << "===========================================" << std::endl
              << "Number of active cells: " << triangulation.n_active_cells()
              << std::endl
              << "Number of degrees of freedom: " << dof_handler.n_dofs()
              << std::endl
              << std::endl;

    constraints.clear();
    DoFTools::make_hanging_node_constraints(dof_handler, constraints);
    constraints.close();

    DynamicSparsityPattern dsp(dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern(dof_handler,
                                    dsp,
                                    constraints,
                                    /*keep_constrained_dofs = */ true);
    sparsity_pattern.copy_from(dsp);

    mass_matrix.reinit(sparsity_pattern);
    matrix.reinit(sparsity_pattern);
    system_matrix.reinit(sparsity_pattern);

    MatrixCreator::create_mass_matrix(dof_handler,
                                      QGauss<dim>(params.fe.quad_order),
                                      mass_matrix);

    solution.reinit(dof_handler.n_dofs());
    system_rhs.reinit(dof_handler.n_dofs());
    rhs.reinit(dof_handler.n_dofs());
    old_rhs.reinit(dof_handler.n_dofs());
  }



  template <int dim>
  void Solver<dim>::assemble_rhs()
  {
    const QGauss<dim> quadrature_formula(params.fe.quad_order);

    rhs = 0;

    FEValues<dim> fe_values(fe,
                            quadrature_formula,
                            update_values | update_quadrature_points |
                              update_JxW_values);

    const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
    const unsigned int n_q_points    = quadrature_formula.size();

    Vector<double> cell_rhs(dofs_per_cell);

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    const double RHS = -material.initialPolymerFraction() *
                       material.polymerRho() * material.dydt(time());

    for (const auto &cell : dof_handler.active_cell_iterators())
      {
        cell_rhs = 0;
        fe_values.reinit(cell);
        for (unsigned int q = 0; q < n_q_points; ++q)
          for (unsigned int i = 0; i < dofs_per_cell; ++i)
            cell_rhs(i) += fe_values.shape_value(i, q) * RHS * fe_values.JxW(q);

        cell->get_dof_indices(local_dof_indices);
        for (unsigned int i = 0; i < dofs_per_cell; ++i)
          rhs(local_dof_indices[i]) += cell_rhs(i);
      }
  }



  template <int dim>
  void Solver<dim>::assemble_matrix()
  {
    const QGauss<dim> quadrature_formula(params.fe.quad_order);

    matrix       = 0;
    max_pressure = 0;

    FEValues<dim> fe_values(fe,
                            quadrature_formula,
                            update_gradients | update_quadrature_points |
                              update_JxW_values);

    const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
    const unsigned int n_q_points    = quadrature_formula.size();

    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    const double K = material.K(time());
    const double D2byP =
      K / material.mu() / material.poreVolumeFraction(time());

    for (const auto &cell : dof_handler.active_cell_iterators())
      {
        cell_matrix = 0;
        fe_values.reinit(cell);
        for (unsigned int q = 0; q < n_q_points; ++q)
          {
            const double D1 = material.D(T);
            const double P  = material.P(time(), solution[q], T);
            const double D2 = D2byP * P;
            if (P > max_pressure)
              max_pressure = P;

            for (unsigned int i = 0; i < dofs_per_cell; ++i)
              for (unsigned int j = 0; j < dofs_per_cell; ++j)
                cell_matrix(i, j) += (D1 + D2) * fe_values.shape_grad(i, q) *
                                     fe_values.shape_grad(j, q) *
                                     fe_values.JxW(q);
          }

        cell->get_dof_indices(local_dof_indices);
        for (unsigned int i = 0; i < dofs_per_cell; ++i)
          for (unsigned int j = 0; j < dofs_per_cell; ++j)
            matrix.add(local_dof_indices[i],
                       local_dof_indices[j],
                       cell_matrix(i, j));
      }
  }



  template <int dim>
  double Solver<dim>::solve_ls()
  {
    SolverControl solver_control(params.ls.max_iter,
                                 params.ls.tol * system_rhs.l2_norm());

    SolverCG<Vector<double>> cg(solver_control);

    PreconditionSSOR<SparseMatrix<double>> preconditioner;
    preconditioner.initialize(system_matrix, 1.0);

    cg.solve(system_matrix, solution, system_rhs, preconditioner);

    constraints.distribute(solution);

    std::cout << "     " << solver_control.last_step()
              << " CG iterations: initial residual = "
              << solver_control.initial_value() / system_rhs.l2_norm()
              << " final residual = "
              << solver_control.last_value() / system_rhs.l2_norm()
              << std::endl;

    return solver_control.initial_value() / system_rhs.l2_norm();
  }



  template <int dim>
  void Solver<dim>::output_results() const
  {
    DataOut<dim> data_out;

    data_out.attach_dof_handler(dof_handler);
    data_out.add_data_vector(solution, "rho");

    data_out.build_patches();

    data_out.set_flags(DataOutBase::VtkFlags(time(), time.step()));

    const std::string filename =
      "solution-" + Utilities::int_to_string(time.step(), 3) + ".vtk";
    std::ofstream output(filename);
    data_out.write_vtk(output);
  }



  template <int dim>
  void Solver<dim>::make_mesh()
  {
    GridGenerator::hyper_cube(triangulation, 0, problem.size);
    triangulation.refine_global(params.mr.j_min);
    const double small = 1e-16 * problem.size;

    for (unsigned int step = 0; step < params.mr.j_max - params.mr.j_min;
         ++step)
      {
        for (auto &cell : triangulation.active_cell_iterators())
          {
            RefinementCase<dim> rf = RefinementCase<dim>::no_refinement;
            for (const auto face_no : cell->face_indices())
              {
                const auto &face = cell->face(face_no);
                if (face->at_boundary())
                  for (int i = 0; i < dim; ++i)
                    if (std::fabs(cell->center()[i] - face->center()[i]) >
                        small)
                      rf = rf | RefinementCase<dim>::cut_axis(i);
              }
            cell->set_refine_flag(rf);
          }

        triangulation.execute_coarsening_and_refinement();
      }
  }



  template <int dim>
  void Solver<dim>::refine_mesh()
  {
    Vector<float> estimated_error_per_cell(triangulation.n_active_cells());

    KellyErrorEstimator<dim>::estimate(
      dof_handler,
      QGauss<dim - 1>(params.fe.quad_order),
      std::map<types::boundary_id, const Function<dim> *>(),
      solution,
      estimated_error_per_cell);

    GridRefinement::refine_and_coarsen_fixed_fraction(triangulation,
                                                      estimated_error_per_cell,
                                                      params.mr.upper,
                                                      params.mr.lower);

    if (triangulation.n_levels() > params.mr.j_max)
      for (const auto &cell :
           triangulation.active_cell_iterators_on_level(params.mr.j_max))
        cell->clear_refine_flag();
    for (const auto &cell :
         triangulation.active_cell_iterators_on_level(params.mr.j_min))
      cell->clear_coarsen_flag();

    std::vector<Vector<double>> transfer_in;
    std::vector<Vector<double>> transfer_out;
    transfer_in.push_back(solution);
    transfer_in.push_back(rhs);

    triangulation.prepare_coarsening_and_refinement();
    SolutionTransfer<dim> solution_trans(dof_handler);
    solution_trans.prepare_for_coarsening_and_refinement(transfer_in);
    triangulation.execute_coarsening_and_refinement();

    setup_system();

    {
      transfer_out.push_back(Vector<double>());
      transfer_out.push_back(Vector<double>());
      transfer_out[0].reinit(dof_handler.n_dofs());
      transfer_out[1].reinit(dof_handler.n_dofs());
    }

    solution_trans.interpolate(transfer_in, transfer_out);
    solution = transfer_out[0];
    rhs      = transfer_out[1];

    constraints.distribute(solution);
    constraints.distribute(rhs);
  }



  template <int dim>
  void Solver<dim>::run()
  {
    make_mesh();
    setup_system();

    // Initial conditions
    VectorTools::interpolate(dof_handler,
                             Functions::ZeroFunction<dim>(),
                             solution);
    output_results();

    assemble_matrix();
    assemble_rhs();

    Vector<double> tmp;
    Vector<double> forcing_terms;

    tmp.reinit(solution.size());
    forcing_terms.reinit(solution.size());

    while (time.loop())
      {
        std::cout << "Time step " << time.step() << " at t=" << time()
                  << std::endl;

        T = problem.T0 + problem.heating_rate * time();

        old_rhs = rhs;
        assemble_rhs();

        unsigned int i = 0;

        // Assemble system_rhs
        forcing_terms = rhs;
        forcing_terms *= time.delta() * time.theta();
        forcing_terms.add(time.delta() * (1 - time.theta()), old_rhs);
        mass_matrix.vmult(system_rhs, solution);
        matrix.vmult(tmp, solution);
        system_rhs.add(-(1 - time.theta()) * time.delta(), tmp);
        system_rhs += forcing_terms;
        constraints.condense(system_rhs);

        // Boundary conditions
        std::map<types::global_dof_index, double> boundary_values;
        VectorTools::interpolate_boundary_values(dof_handler,
                                                 0,
                                                 Functions::ZeroFunction<dim>(),
                                                 boundary_values);

        do
          {
            assemble_matrix();

            // Assemble system_matrix
            system_matrix.copy_from(mass_matrix);
            system_matrix.add(time.theta() * time.delta(), matrix);
            constraints.condense(system_matrix);

            MatrixTools::apply_boundary_values(boundary_values,
                                               system_matrix,
                                               solution,
                                               system_rhs);
          }
        while (solve_ls() > params.ns.tol && ++i < params.ns.max_iter);

        output_results();
        std::cout << "y = " << material.y(time()) << " T = " << T
                  << " max(p) = " << max_pressure << std::endl;

        if (time.step() % params.mr.n_steps == 0)
          {
            refine_mesh();
            tmp.reinit(solution.size());
            forcing_terms.reinit(solution.size());
          }
      }
  }
} // namespace ThermalDebinding


int main()
{
  try
    {
      using namespace ThermalDebinding;

      Parameters parameters("parameters.prm");
      Solver<2>  solver(parameters);
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
