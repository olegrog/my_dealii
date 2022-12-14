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

#include "Parameters.h"

namespace ThermalDebinding
{
  Problem::Problem()
    : ParameterAcceptor("Problem")
  {
    add_parameter("Hyper cube size", size = 1e-2);
    add_parameter("Initial temperature", T0 = 300);
    add_parameter("Heating rate", heating_rate = 1.667e-3);
    add_parameter("Minimum pressure", min_pressure = 1e3);
  }



  FiniteElements::FiniteElements()
    : ParameterAcceptor("Finite elements")
  {
    add_parameter("Polynomial degree", poly_degree = 1);
    add_parameter("Renumbering", renumbering = false);
    add_parameter("Lumped Mass Matrix", lumped_mass_matrix = false);
  }



  MeshRefinement::MeshRefinement()
    : ParameterAcceptor("Mesh refinement")
  {
    add_parameter("Minimum level", j_min = 4);
    add_parameter("Maximum level", j_max = 8);
    add_parameter("Number of time steps", n_steps = 1);
    add_parameter("Refining threshold", upper = 0.5);
    add_parameter("Coarsening threshold", lower = 0);
    add_parameter("Maximum number of cells", max_n_cells = 5000);
  }



  LinearSolver::LinearSolver()
    : ParameterAcceptor("Linear solver")
  {
    add_parameter("Max iterations", max_iter = 100);
    add_parameter("Tolerance", tol = 1e-9);
    add_parameter("Reduce", reduce = 1e-2);
    add_parameter("Preconditioner relaxation", preconditioner_relax = 1.0);
  }



  NonlinearSolver::NonlinearSolver()
    : ParameterAcceptor("Nonlinear solver")
  {
    add_parameter("Max iterations", max_iter = 15);
    add_parameter("Tolerance", tol = 1e-8);
  }



  Output::Output()
    : ParameterAcceptor("Output")
  {
    add_parameter("Write VTK files", write_vtk_files = true);
    add_parameter("Write mesh", write_mesh = false);
    add_parameter("Mesh format",
                  mesh_format_str = "vtk",
                  "",
                  ParameterAcceptor::prm,
                  Patterns::Selection(GridOut::get_output_format_names()));
    add_parameter("Number of time steps", n_steps = 1);
    add_parameter("Log verbosity", verbosity = 0);
    add_parameter("Profiling", profiling = false);
  }



  void Output::parse_parameters(ParameterHandler &)
  {
    mesh_format = GridOut::parse_output_format(mesh_format_str);
  }

} // namespace ThermalDebinding
