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
#include <deal.II/base/timer.h>

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
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_q.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/solution_transfer.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/derivative_approximation.h>

#include <fstream>
#include <iostream>
#include <algorithm>

#include "Parameters.h"
#include "Time.h"

namespace FractureHealing
{
  using namespace dealii;


  template <int dim>
  class Solver
  {
  public:
    Solver(Parameters<dim> &parameters, Time &time);
    void run();

  private:
    void   setup_system();
    void   assemble_rhs();
    void   assemble_matrix();
    double solve_ls();
    void   output_results();
    void   make_mesh();
    bool   refine_mesh();
    void   transfer_solution();

    Parameters<dim> &params;
    Time &           time;

    Triangulation<dim> triangulation;
    DoFHandler<dim>    dof_handler;
    FESystem<dim>      fe;

    AffineConstraints<double> constraints;

    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> mass_matrix;
    SparseMatrix<double> matrix;
    SparseMatrix<double> system_matrix;

    Vector<double> solution;
    Vector<double> old_rhs;
    Vector<double> rhs;
    Vector<double> system_rhs;

    TimerOutput computing_timer;
  };



  template <int dim>
  Solver<dim>::Solver(Parameters<dim> &parameters, Time &time)
    : params(parameters)
    , time(time)
    , dof_handler(triangulation)
    , fe(FE_Q<dim>(params.fe.poly_degree), Model::n_components)
    , computing_timer(std::cout,
                      TimerOutput::never,
                      TimerOutput::cpu_and_wall_times_grouped)
  {}



  template <int dim>
  void Solver<dim>::setup_system()
  {
    if (params.output.write_mesh)
      {
        Timer timer;
        if (params.output.verbosity > 0)
          std::cout << " -- Writing the mesh... " << std::flush;
        TimerOutput::Scope t(computing_timer, "Writing the mesh");

        std::ofstream out("grid-" + Utilities::int_to_string(time.step(), 3) +
                          GridOut::default_suffix(params.output.mesh_format));
        GridOut       grid_out;
        grid_out.write(triangulation, out, params.output.mesh_format);

        if (params.output.verbosity > 0)
          std::cout << "done (" << timer.cpu_time() << "s)" << std::endl;
      }

    dof_handler.distribute_dofs(fe);

    if (params.fe.renumbering)
      {
        Timer timer;
        if (params.output.verbosity > 0)
          std::cout << " -- Renumbering DoFs... " << std::flush;
        TimerOutput::Scope t(computing_timer, "Renumbering DoFs");

        DoFRenumbering::Cuthill_McKee(dof_handler);

        if (params.output.verbosity > 0)
          std::cout << "done (" << timer.cpu_time() << "s)" << std::endl;
      }

    std::cout << std::endl
              << "===========================================" << std::endl
              << "Number of active cells: " << triangulation.n_active_cells()
              << std::endl
              << "Number of degrees of freedom: " << dof_handler.n_dofs()
              << std::endl;

    Timer timer;
    if (params.output.verbosity > 0)
      std::cout << " -- Setting up system... " << std::flush;
    TimerOutput::Scope t(computing_timer, "Setting up system");

    {
      std::map<types::boundary_id, unsigned int> id_faces;

      // NB: loop over all the cells, not only active
      for (auto &cell : triangulation.cell_iterators())
        for (auto &face : cell->face_iterators())
          if (face->at_boundary())
            {
              const types::boundary_id id = params.bc.get_id(face->center());
              id_faces[id]++;
              face->set_boundary_id(id);
            }

      for (auto [key, value] : id_faces)
        if (key < Model::n_components)
          std::cout << "Boundary #" << key << " contains " << value << " faces"
                    << std::endl;
        else
          std::cout << "Number of faces without specified BC: " << value
                    << std::endl;

      std::cout << std::endl;
    }

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
                                      QGauss<dim>(fe.degree + 1),
                                      mass_matrix);

    solution.reinit(dof_handler.n_dofs());
    system_rhs.reinit(dof_handler.n_dofs());
    rhs.reinit(dof_handler.n_dofs());
    old_rhs.reinit(dof_handler.n_dofs());

    if (params.output.verbosity > 0)
      {
        std::cout << "done (" << timer.cpu_time() << "s)" << std::endl;
        unsigned int n_nonzero = sparsity_pattern.n_nonzero_elements();
        std::cout << "Sparsity pattern: n_nonzero = " << n_nonzero << " = "
                  << std::setprecision(3)
                  << 100. * n_nonzero / std::pow(dof_handler.n_dofs(), 2) << "%"
                  << std::endl;
      }

    if (params.output.verbosity > 1)
      {
        std::ofstream out("sparsity_pattern.svg");
        sparsity_pattern.print_svg(out);
      }
  }



  template <int dim>
  void Solver<dim>::assemble_rhs()
  {
    Timer timer;
    if (params.output.verbosity > 0)
      std::cout << " -- Assembling the RHS... " << std::flush;
    TimerOutput::Scope t(computing_timer, "Assembling the RHS");

    const QGauss<dim> quadrature_formula(fe.degree + 1);

    rhs = 0;

    FEValues<dim> fe_values(fe,
                            quadrature_formula,
                            update_values | update_JxW_values);

    const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
    const unsigned int n_q_points    = quadrature_formula.size();

    Vector<double> cell_rhs(dofs_per_cell);

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    std::vector<Model::vector_type> Q_values(n_q_points);
    std::vector<Vector<double>>     solution_values(
      n_q_points, Vector<double>(Model::n_components));

    for (const auto &cell : dof_handler.active_cell_iterators())
      {
        cell_rhs = 0;
        fe_values.reinit(cell);
        fe_values.get_function_values(solution, solution_values);
        params.model.get_source_vectors(Q_values, solution_values);

        for (const unsigned int i : fe_values.dof_indices())
          {
            const unsigned int component_i =
              fe.system_to_component_index(i).first;

            for (const unsigned int q_point :
                 fe_values.quadrature_point_indices())
              cell_rhs(i) += fe_values.shape_value(i, q_point) *
                             Q_values[q_point][component_i] *
                             fe_values.JxW(q_point);
          }

        cell->get_dof_indices(local_dof_indices);
        for (const unsigned int i : fe_values.dof_indices())
          rhs(local_dof_indices[i]) += cell_rhs(i);
      }

    if (params.output.verbosity > 0)
      std::cout << "done (" << timer.cpu_time() << "s)" << std::endl;
  }



  template <int dim>
  void Solver<dim>::assemble_matrix()
  {
    Timer timer;
    if (params.output.verbosity > 0)
      std::cout << " -- Assembling the matrix... " << std::flush;
    TimerOutput::Scope t(computing_timer, "Assembling the matrix");

    const QGauss<dim> quadrature_formula(fe.degree + 1);

    matrix = 0;

    FEValues<dim> fe_values(fe,
                            quadrature_formula,
                            update_values | update_gradients |
                              update_JxW_values);

    const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
    const unsigned int n_q_points    = quadrature_formula.size();

    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    std::vector<Model::matrix_type> D_values(n_q_points);
    std::vector<Model::matrix_type> R_values(n_q_points);
    std::vector<Vector<double>>     solution_values(
      n_q_points, Vector<double>(Model::n_components));

    std::vector<double>         phi(dofs_per_cell);
    std::vector<Tensor<1, dim>> grad_phi(dofs_per_cell);

    for (const auto &cell : dof_handler.active_cell_iterators())
      {
        cell_matrix = 0;
        fe_values.reinit(cell);
        fe_values.get_function_values(solution, solution_values);
        params.model.get_diffusion_matrices(D_values, solution_values);
        params.model.get_reaction_matrices(R_values, solution_values);

        for (const unsigned int q_point : fe_values.quadrature_point_indices())
          {
            for (const unsigned int k : fe_values.dof_indices())
              {
                phi[k]      = fe_values.shape_value(k, q_point);
                grad_phi[k] = fe_values.shape_grad(k, q_point);
              }

            for (const unsigned int i : fe_values.dof_indices())
              {
                const unsigned int component_i =
                  fe.system_to_component_index(i).first;

                for (const unsigned int j : fe_values.dof_indices())
                  {
                    const unsigned int component_j =
                      fe.system_to_component_index(j).first;

                    cell_matrix(i, j) +=
                      (D_values[q_point][component_i][component_j] *
                         grad_phi[i] * grad_phi[j] -
                       R_values[q_point][component_i][component_j] * phi[i] *
                         phi[j]) *
                      fe_values.JxW(q_point);
                  }
              }
          }

        cell->get_dof_indices(local_dof_indices);
        for (const unsigned int i : fe_values.dof_indices())
          for (const unsigned int j : fe_values.dof_indices())
            matrix.add(local_dof_indices[i],
                       local_dof_indices[j],
                       cell_matrix(i, j));
      }

    if (params.output.verbosity > 0)
      std::cout << "done (" << timer.cpu_time() << "s)" << std::endl;
  }



  template <int dim>
  double Solver<dim>::solve_ls()
  {
    Timer timer;
    if (params.output.verbosity > 0)
      std::cout << " -- Solving a linear system... " << std::flush;
    TimerOutput::Scope t(computing_timer, "Solving a linear system");

    ReductionControl solver_control(params.ls.max_iter,
                                    params.ls.tol * system_rhs.l2_norm(),
                                    params.ls.reduce);

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

    return solver_control.initial_value() / system_rhs.l2_norm();
  }



  template <int dim>
  void Solver<dim>::output_results()
  {
    if (!params.output.write_vtk_files)
      return;

    if (time.step() % params.output.n_steps != 0)
      return;

    Timer timer;
    if (params.output.verbosity > 0)
      std::cout << " -- Writing results... " << std::flush;
    TimerOutput::Scope t(computing_timer, "Writing results");

    DataOut<dim> data_out;

    data_out.attach_dof_handler(dof_handler);
    {
      // Remove all spaces from Model::component_names
      std::vector<std::string> component_names = Model::component_names;
      for (auto &name : component_names)
        std::replace(name.begin(), name.end(), ' ', '_');
      data_out.add_data_vector(solution, component_names);
    }
    data_out.build_patches();

    DataOutBase::VtkFlags output_flags(time(), time.step());
    output_flags.compression_level =
      DataOutBase::VtkFlags::ZlibCompressionLevel::best_speed;
    data_out.set_flags(output_flags);

    const std::string filename =
      "solution-" + Utilities::int_to_string(time.step(), 3) + ".vtu";
    std::ofstream output(filename);
    data_out.write_vtu(output);

    if (params.output.verbosity > 0)
      std::cout << "done (" << timer.cpu_time() << "s)" << std::endl;
  }



  template <int dim>
  void Solver<dim>::make_mesh()
  {
    Timer timer;
    if (params.output.verbosity > 0)
      std::cout << " -- Making a mesh... " << std::flush;
    TimerOutput::Scope t(computing_timer, "Making a mesh");

    const Geometry<dim> &geo     = params.geometry;
    const unsigned int   n_steps = params.mr.j_max - params.mr.j_min;

    GridGenerator::subdivided_hyper_L(triangulation,
                                      geo.repetitions,
                                      Point<dim>(),
                                      geo.size,
                                      geo.n_cells_to_remove);
    triangulation.refine_global(params.mr.j_min);

    // Adapt mesh to the initial values
    for (unsigned int step = 0; step < n_steps; ++step)
      {
        std::cout << "Number of active cells: "
                  << triangulation.n_active_cells() << std::endl;

        dof_handler.distribute_dofs(fe);
        solution.reinit(dof_handler.n_dofs());
        VectorTools::interpolate(dof_handler, params.ic, solution);

        Vector<float> gradient_indicator(triangulation.n_active_cells());

        DerivativeApproximation::approximate_gradient(dof_handler,
                                                      solution,
                                                      gradient_indicator);

        unsigned int cell_no = 0;
        for (const auto &cell : dof_handler.active_cell_iterators())
          gradient_indicator(cell_no++) *=
            std::pow(cell->diameter(), 1 + dim / 2.);

        GridRefinement::refine_and_coarsen_fixed_fraction(
          triangulation,
          gradient_indicator,
          params.mr.ic_upper,
          0,
          params.mr.max_n_cells);

        triangulation.execute_coarsening_and_refinement();
      }

    if (params.output.verbosity > 0)
      std::cout << "done (" << timer.cpu_time() << "s)" << std::endl;
  }



  template <int dim>
  bool Solver<dim>::refine_mesh()
  {
    Timer timer;
    if (params.output.verbosity > 0)
      std::cout << " -- Refining the mesh... " << std::flush;
    TimerOutput::Scope t(computing_timer, "Refining the mesh");

    Vector<float> estimated_error_per_cell(triangulation.n_active_cells());

    KellyErrorEstimator<dim>::estimate(
      dof_handler,
      QGauss<dim - 1>(fe.degree + 1),
      std::map<types::boundary_id, const Function<dim> *>(),
      solution,
      estimated_error_per_cell);

    GridRefinement::refine_and_coarsen_fixed_fraction(triangulation,
                                                      estimated_error_per_cell,
                                                      params.mr.upper,
                                                      params.mr.lower,
                                                      params.mr.max_n_cells);

    if (triangulation.n_levels() > params.mr.j_max)
      for (const auto &cell :
           triangulation.active_cell_iterators_on_level(params.mr.j_max))
        cell->clear_refine_flag();
    for (const auto &cell :
         triangulation.active_cell_iterators_on_level(params.mr.j_min))
      cell->clear_coarsen_flag();

    bool is_changed = triangulation.prepare_coarsening_and_refinement();

    if (params.output.verbosity > 0)
      std::cout << "done (" << timer.cpu_time() << "s)" << std::endl;

    return is_changed;
  }



  template <int dim>
  void Solver<dim>::transfer_solution()
  {
    const Vector<double>  previous_solution = solution;
    SolutionTransfer<dim> solution_trans(dof_handler);

    solution_trans.prepare_for_coarsening_and_refinement(previous_solution);
    triangulation.execute_coarsening_and_refinement();

    setup_system();

    solution_trans.interpolate(previous_solution, solution);
    constraints.distribute(solution);
  }



  template <int dim>
  void Solver<dim>::run()
  {
    make_mesh();
    setup_system();

    // Initial conditions
    if (params.ic.project_functions)
      VectorTools::project(dof_handler,
                           constraints,
                           QGauss<dim>(fe.degree + 1),
                           params.ic,
                           solution);
    else
      VectorTools::interpolate(dof_handler, params.ic, solution);

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

        computing_timer.reset();

        old_rhs = rhs;
        assemble_rhs();

        {
          Timer timer;
          if (params.output.verbosity > 0)
            std::cout << " -- Assembling the system RHS... " << std::flush;
          TimerOutput::Scope t(computing_timer, "Assembling the system RHS");

          forcing_terms = rhs;
          forcing_terms *= time.delta() * time.theta();
          forcing_terms.add(time.delta() * (1 - time.theta()), old_rhs);
          mass_matrix.vmult(system_rhs, solution);
          matrix.vmult(tmp, solution);
          system_rhs.add(-(1 - time.theta()) * time.delta(), tmp);
          system_rhs += forcing_terms;
          constraints.condense(system_rhs);

          if (params.output.verbosity > 0)
            std::cout << "done (" << timer.cpu_time() << "s)" << std::endl;
        }

        // Boundary conditions
        std::map<types::global_dof_index, double> boundary_values;
        params.bc.interpolate_boundary_values(dof_handler, boundary_values);

        unsigned int i = 0;

        do
          {
            assemble_matrix();

            Timer timer;
            if (params.output.verbosity > 0)
              std::cout << " -- Assembling the system matrix... " << std::flush;
            TimerOutput::Scope t(computing_timer,
                                 "Assembling the system matrix");

            system_matrix.copy_from(mass_matrix);
            system_matrix.add(time.theta() * time.delta(), matrix);
            constraints.condense(system_matrix);

            // NB: last `false` makes matrix unsymmetric, but CG still works!
            MatrixTools::apply_boundary_values(boundary_values,
                                               system_matrix,
                                               solution,
                                               system_rhs,
                                               /* eliminate_columns = */ false);

            if (params.output.verbosity > 0)
              std::cout << "done (" << timer.cpu_time() << "s)" << std::endl;
          }
        while (solve_ls() > params.ns.tol && ++i < params.ns.max_iter);

        constraints.distribute(solution);

        output_results();

        if (time.step() % params.mr.n_steps == 0)
          {
            if (refine_mesh())
              {
                transfer_solution();
                tmp.reinit(solution.size());
                forcing_terms.reinit(solution.size());
                assemble_rhs();
              }
          }

        if (params.output.profiling)
          computing_timer.print_summary();
      }
  }
} // namespace FractureHealing



int main(int argc, char *argv[])
{
  try
    {
      using namespace FractureHealing;

      constexpr int dim = 2;

      Parameters<dim> parameters;
      Time            time;

      const std::string prm_file = argc > 1 ? argv[1] : "parameters.prm";
      std::cout << "Reading parameters... " << std::flush;
      ParameterAcceptor::initialize(prm_file);
      std::cout << "done" << std::endl;

      deallog.depth_console(parameters.output.verbosity);

      Solver<dim> solver(parameters, time);
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
