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

#include <deal.II/base/point.h>
#include <deal.II/base/parameter_acceptor.h>
#include <deal.II/grid/grid_out.h>

namespace BimetallicStrip
{
  using namespace dealii;

  template <int dim>
  struct Problem : ParameterAcceptor
  {
    Problem();

    void parse_parameters(ParameterHandler &prm) override;

    Point<dim>   dimensions;
    double       alpha1, alpha2;
    double       E, nu;
    double       mu, lambda, threeK;
    unsigned int axis;
  };

  struct FiniteElements : ParameterAcceptor
  {
    FiniteElements();

    unsigned int poly_degree;
  };

  struct MeshRefinement : ParameterAcceptor
  {
    MeshRefinement();

    unsigned int j_min;
    unsigned int j_max;
    double       upper;
    double       lower;
    double       max_n_cells;
  };

  struct LinearSolver : ParameterAcceptor
  {
    LinearSolver();

    unsigned int max_iter;
    double       tol;
    double       preconditioner_relax;
  };

  struct Output : ParameterAcceptor
  {
    Output();

    unsigned int verbosity;
  };

  template <int dim>
  struct Parameters
  {
    Problem<dim>   problem;
    FiniteElements fe;
    MeshRefinement mr;
    LinearSolver   ls;
    Output         output;
  };

} // namespace BimetallicStrip
