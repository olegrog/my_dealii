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
#include <deal.II/dofs/dof_renumbering.h>
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
#include "Material.h"
#include "Time.h"
#include "PostProcess.h"

namespace ThermalDebinding
{
  using namespace dealii;


  template <int dim>
  class Solver
  {
  public:
    Solver(Parameters &parameters, Material &material, Time &time);
    void run();

  private:
    void   setup_system();
    void   assemble_rhs();
    void   assemble_matrix();
    void   find_extreme_values() const;
    double solve_ls();
    void   output_results() const;
    void   make_mesh();
    bool   refine_mesh();

    Parameters &   params;
    const Problem &problem;
    Material &     material;
    Time &         time;

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

    double T; // current temperature [K]
  };



  template <int dim>
  Solver<dim>::Solver(Parameters &parameters, Material &material, Time &time)
    : params(parameters)
    , problem(params.problem)
    , material(material)
    , time(time)
    , fe(params.fe.poly_degree)
    , dof_handler(triangulation)
    , T(problem.T0)
  {}



  template <int dim>
  void Solver<dim>::setup_system()
  {
    if (params.output.verbosity > 0)
      std::cout << " -- Setup a system" << std::endl;

    if (params.output.write_mesh)
      {
        std::cout << "Writing the mesh... " << std::flush;
        std::ofstream out("grid-" + Utilities::int_to_string(time.step(), 3) +
                          GridOut::default_suffix(params.output.mesh_format));
        GridOut       grid_out;
        grid_out.write(triangulation, out, params.output.mesh_format);
        std::cout << "done" << std::endl;
      }

    dof_handler.distribute_dofs(fe);

    if (params.fe.renumbering)
      {
        std::cout << "Renumbering DoFs... " << std::flush;
        DoFRenumbering::Cuthill_McKee(dof_handler);
        std::cout << "done" << std::endl;
      }

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

    if (params.output.verbosity > 1)
      {
        unsigned int n_nonzero = sparsity_pattern.n_nonzero_elements();
        std::cout << "Sparsity pattern: n_nonzero = " << n_nonzero << " = "
                  << std::setprecision(3)
                  << 100. * n_nonzero / std::pow(dof_handler.n_dofs(), 2)
                  << "%\n";
        std::ofstream out("sparsity_pattern.svg");
        sparsity_pattern.print_svg(out);
      }

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
    if (params.output.verbosity > 0)
      std::cout << " -- Assemble the RHS" << std::endl;

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
                       material.polymerRho() * material.dydt(T);

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
    if (params.output.verbosity > 0)
      std::cout << " -- Assemble the matrix" << std::endl;

    const QGauss<dim> quadrature_formula(params.fe.quad_order);

    matrix = 0;

    FEValues<dim> fe_values(fe,
                            quadrature_formula,
                            update_values | update_gradients |
                              update_quadrature_points | update_JxW_values);

    const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
    const unsigned int n_q_points    = quadrature_formula.size();

    FullMatrix<double>  cell_matrix(dofs_per_cell, dofs_per_cell);
    std::vector<double> solution_values(n_q_points);

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    std::vector<Tensor<1, dim>> grad_phi(dofs_per_cell);

    for (const auto &cell : dof_handler.active_cell_iterators())
      {
        cell_matrix = 0;
        fe_values.reinit(cell);
        fe_values.get_function_values(solution, solution_values);

        for (unsigned int q = 0; q < n_q_points; ++q)
          {
            const double D1 = material.D1(T);
            const double P  = material.P(solution_values[q], T);
            const double D2 = material.D2(P);

            for (unsigned int k = 0; k < dofs_per_cell; ++k)
              grad_phi[k] = fe_values.shape_grad(k, q);

            for (unsigned int i = 0; i < dofs_per_cell; ++i)
              for (unsigned int j = 0; j <= i; ++j)
                cell_matrix(i, j) +=
                  (D1 + D2) * grad_phi[i] * grad_phi[j] * fe_values.JxW(q);
          }

        cell->get_dof_indices(local_dof_indices);
        for (unsigned int i = 0; i < dofs_per_cell; ++i)
          for (unsigned int j = 0; j < dofs_per_cell; ++j)
            if (j <= i)
              matrix.add(local_dof_indices[i],
                         local_dof_indices[j],
                         cell_matrix(i, j));
            else
              matrix.add(local_dof_indices[i],
                         local_dof_indices[j],
                         cell_matrix(j, i));
      }
  }



  template <int dim>
  void Solver<dim>::find_extreme_values() const
  {
    if (params.output.verbosity > 0)
      std::cout << " -- Find the extreme values" << std::endl;

    const QIterated<dim> quadrature_formula(QTrapezoid<1>(),
                                            params.fe.quad_order);
    FEValues<dim>        fe_values(fe, quadrature_formula, update_values);

    const unsigned int n_q_points = quadrature_formula.size();

    std::vector<double> solution_values(n_q_points);

    double maxP   = 0;
    double maxRho = 0;
    double minD1  = std::numeric_limits<double>::max();
    double maxD1  = 0;
    double minD2  = std::numeric_limits<double>::max();
    double maxD2  = 0;

    for (const auto &cell : dof_handler.active_cell_iterators())
      {
        fe_values.reinit(cell);
        fe_values.get_function_values(solution, solution_values);

        for (unsigned int q = 0; q < n_q_points; ++q)
          {
            const double P  = material.P(solution_values[q], T);
            const double D1 = material.D1(T);
            const double D2 = material.D2(P);

            maxP   = std::max(maxP, P);
            maxRho = std::max(maxRho, solution_values[q]);
            minD1  = std::min(minD1, D1);
            maxD1  = std::max(maxD1, D1);
            minD2  = std::min(minD2, D2);
            maxD2  = std::max(maxD2, D2);
          }
      }

    std::cout << "max(P) = " << maxP << " max(Rho) = " << maxRho
              << " min(D1) = " << minD1 << " max(D1) = " << maxD1
              << " min(D2) = " << minD2 << " max(D2) = " << maxD2 << std::endl;
  }



  template <int dim>
  double Solver<dim>::solve_ls()
  {
    if (params.output.verbosity > 0)
      std::cout << " -- Solve a linear system" << std::endl;

    ReductionControl solver_control(params.ls.max_iter,
                                    params.ls.tol * system_rhs.l2_norm(),
                                    params.ls.reduce);

    SolverCG<Vector<double>> cg(solver_control);

    PreconditionSSOR<SparseMatrix<double>> preconditioner;
    preconditioner.initialize(system_matrix, params.ls.preconditioner_relax);

    cg.solve(system_matrix, solution, system_rhs, preconditioner);

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
    if (params.output.verbosity > 0)
      std::cout << " -- Output the results" << std::endl;

    if (!params.output.write_vtk_files)
      return;

    if (time.step() % params.output.n_steps == 0)
      {
        DataOut<dim>         data_out;
        ComputePressure<dim> compute_pressure(material, T);

        data_out.attach_dof_handler(dof_handler);
        data_out.add_data_vector(solution, "density");
        data_out.add_data_vector(solution, compute_pressure);

        data_out.build_patches();

        data_out.set_flags(DataOutBase::VtkFlags(time(), time.step()));

        const std::string filename =
          "solution-" + Utilities::int_to_string(time.step(), 3) + ".vtk";
        std::ofstream output(filename);
        data_out.write_vtk(output);
      }
  }



  template <int dim>
  void Solver<dim>::make_mesh()
  {
    if (params.output.verbosity > 0)
      std::cout << " -- Make a mesh" << std::endl;

    GridGenerator::hyper_cube(triangulation, 0, problem.size);
    triangulation.refine_global(params.mr.j_min);
    const double small = 1e-15 * problem.size;

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
  bool Solver<dim>::refine_mesh()
  {
    // Anisotropic refinement does not work in the 3D case
    if (dim == 3)
      return false;

    if (params.output.verbosity > 0)
      std::cout << " -- Refine the mesh" << std::endl;

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

    if (is_changed)
      {
        Vector<double>        previous_solution = solution;
        SolutionTransfer<dim> solution_trans(dof_handler);
        solution_trans.prepare_for_coarsening_and_refinement(previous_solution);
        triangulation.execute_coarsening_and_refinement();
        setup_system();
        solution_trans.interpolate(previous_solution, solution);
        constraints.distribute(solution);
      }

    return is_changed;
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
        material.evolve(T, time.delta());

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

        constraints.distribute(solution);

        output_results();

        std::cout << "T = " << T
                  << " porosity = " << material.poreVolumeFraction();
        for (unsigned int i = 0; i < material.species().size(); i++)
          std::cout << " y" << i + 1 << " = " << material.species()[i].y;
        std::cout << std::endl;

        find_extreme_values();

        if (time.step() % params.mr.n_steps == 0)
          {
            if (refine_mesh())
              {
                tmp.reinit(solution.size());
                forcing_terms.reinit(solution.size());
                assemble_rhs();
              }
          }
      }
  }
} // namespace ThermalDebinding



int main(int argc, char *argv[])
{
  try
    {
      using namespace ThermalDebinding;

      Parameters parameters;
      Material   material;
      Time       time;

      const std::string prm_file = argc > 1 ? argv[1] : "parameters.prm";
      std::cout << "Reading parameters... " << std::flush;
      ParameterAcceptor::initialize(prm_file);
      std::cout << "done" << std::endl;

      deallog.depth_console(parameters.output.verbosity);

      Solver<3> solver(parameters, material, time);
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
