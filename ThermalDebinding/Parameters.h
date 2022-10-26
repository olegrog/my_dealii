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

#include <deal.II/base/parameter_acceptor.h>
#include <deal.II/grid/grid_out.h>

namespace ThermalDebinding
{
  using namespace dealii;

  struct Problem : ParameterAcceptor
  {
    Problem();

    double size;
    double T0;
    double heating_rate;
  };

  struct FiniteElements : ParameterAcceptor
  {
    FiniteElements();

    unsigned int poly_degree;
    bool         renumbering;
    bool         lumped_mass_matrix;
  };

  struct MeshRefinement : ParameterAcceptor
  {
    MeshRefinement();

    unsigned int j_min;
    unsigned int j_max;
    unsigned int n_steps;
    double       upper;
    double       lower;
    double       max_n_cells;
  };

  struct LinearSolver : ParameterAcceptor
  {
    LinearSolver();

    unsigned int max_iter;
    double       tol;
    double       reduce;
    double       preconditioner_relax;
  };

  struct NonlinearSolver : ParameterAcceptor
  {
    NonlinearSolver();

    unsigned int max_iter;
    double       tol;
  };

  struct Output : ParameterAcceptor
  {
    Output();

    void parse_parameters(ParameterHandler &prm) override;

    bool                  write_vtk_files;
    bool                  write_mesh;
    std::string           mesh_format_str;
    GridOut::OutputFormat mesh_format;
    unsigned int          n_steps;
    unsigned int          verbosity;
  };

  struct Parameters
  {
    Problem         problem;
    FiniteElements  fe;
    MeshRefinement  mr;
    LinearSolver    ls;
    NonlinearSolver ns;
    Output          output;
  };

} // namespace ThermalDebinding
