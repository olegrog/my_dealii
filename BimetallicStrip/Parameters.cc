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

namespace BimetallicStrip
{
  template <int dim>
  Problem<dim>::Problem()
    : ParameterAcceptor("Problem")
  {
    add_parameter("Dimensions", dimensions);
    add_parameter("Young's modulus", E = 1);
    add_parameter("Poisson's ratio", nu = 0.32);
    add_parameter("Alpha1", alpha1 = 0.025);
    add_parameter("Alpha2", alpha2 = 0);
    add_parameter("Axis", axis = 0);
  }



  template <int dim>
  void Problem<dim>::parse_parameters(ParameterHandler &)
  {
    mu     = E / (2 * (1 + nu));
    lambda = nu * E / (1 + nu) / (1 - 2 * nu);
    threeK = E / (1 - 2 * nu);
  }



  // Explicit instantiation
  template class Problem<1>;
  template class Problem<2>;
  template class Problem<3>;



  FiniteElements::FiniteElements()
    : ParameterAcceptor("Finite elements")
  {
    add_parameter("Polynomial degree", poly_degree = 1);
  }



  MeshRefinement::MeshRefinement()
    : ParameterAcceptor("Mesh refinement")
  {
    add_parameter("Minimum level", j_min = 4);
    add_parameter("Maximum level", j_max = 8);
    add_parameter("Refining threshold", upper = 0.5);
    add_parameter("Coarsening threshold", lower = 0);
    add_parameter("Maximum number of cells", max_n_cells = 5000);
  }



  LinearSolver::LinearSolver()
    : ParameterAcceptor("Linear solver")
  {
    add_parameter("Max iterations", max_iter = 100);
    add_parameter("Tolerance", tol = 1e-8);
    add_parameter("Preconditioner relaxation", preconditioner_relax = 1.0);
  }



  Output::Output()
    : ParameterAcceptor("Output")
  {
    add_parameter("Log verbosity", verbosity = 0);
  }

} // namespace BimetallicStrip
