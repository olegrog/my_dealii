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
#include <deal.II/base/function_parser.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/dofs/dof_handler.h>

namespace FractureHealing
{
  using namespace dealii;

  struct Model : ParameterAcceptor
  {
    static constexpr unsigned int n_components = 7;

    static const std::vector<std::string> component_names;

    using vector_type = Tensor<1, n_components>;
    using matrix_type = Tensor<2, n_components>;

    Model();

    void
    get_diffusion_matrices(std::vector<matrix_type> &         matrices,
                           const std::vector<Vector<double>> &values) const;

    void get_reaction_matrices(std::vector<matrix_type> &         matrices,
                               const std::vector<Vector<double>> &values) const;

    void get_source_vectors(std::vector<vector_type> &         vectors,
                            const std::vector<Vector<double>> &values) const;

    double Dh, Kh, Ck, Kk, Am0, Km, Ac0, Kc, Ab0, Kb, alpha_m, alpha_c, alpha_b,
      Y1, H1, Y2, H2, Y3, H3, Bec, db, Pcs, Pbs, Qcd1, Qcd2, Qbd, Dgc, Dgb, Ggc,
      Hgc, Kgc, Ggb, Hgb, dgc, dgb;
  };

  template <int dim>
  struct Geometry : ParameterAcceptor
  {
    Geometry();

    std::vector<unsigned int> repetitions;
    Point<dim>                size;
    std::vector<int>          n_cells_to_remove;
  };

  template <int dim>
  struct InitialValues : ParameterAcceptor, FunctionParser<dim>
  {
    InitialValues();

    void declare_parameters(ParameterHandler &prm) override;
    void parse_parameters(ParameterHandler &prm) override;

    bool project_functions;
  };

  template <int dim>
  class BoundaryValues : ParameterAcceptor
  {
  public:
    BoundaryValues();

    void declare_parameters(ParameterHandler &prm) override;
    void parse_parameters(ParameterHandler &prm) override;

    types::boundary_id get_id(const Point<dim> &point) const;

    void interpolate_boundary_values(
      const DoFHandler<dim> &                    dof_handler,
      std::map<types::global_dof_index, double> &boundary_values) const;

  private:
    struct Boundary
    {
      Boundary(const std::string &             location_expression,
               const std::vector<std::string> &value_expressions,
               const std::vector<bool> &       cmask);

      const FunctionParser<dim> location;
      FunctionParser<dim>       values;
      const ComponentMask       component_mask;
    };

    static constexpr unsigned int          max_n_boundaries = 10;
    std::map<types::boundary_id, Boundary> boundaries;

    std::string subsection(types::boundary_id id) const
    {
      return "Boundary_" + Utilities::int_to_string(id);
    }
  };

  struct FiniteElements : ParameterAcceptor
  {
    FiniteElements();

    unsigned int poly_degree;
    bool         renumbering;
  };

  struct MeshRefinement : ParameterAcceptor
  {
    MeshRefinement();

    unsigned int j_min;
    unsigned int j_max;
    unsigned int n_steps;
    double       ic_upper;
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
    bool                  profiling;
  };

  template <int dim>
  struct Parameters
  {
    Model               model;
    Geometry<dim>       geometry;
    InitialValues<dim>  ic;
    BoundaryValues<dim> bc;
    FiniteElements      fe;
    MeshRefinement      mr;
    LinearSolver        ls;
    NonlinearSolver     ns;
    Output              output;
  };

  template <int dim>
  inline types::boundary_id
  BoundaryValues<dim>::get_id(const Point<dim> &point) const
  {
    types::boundary_id id = max_n_boundaries;
    for (const auto &[boundary_id, boundary] : boundaries)
      if (boundary.location.value(point))
        {
          Assert(id == max_n_boundaries,
                 ExcMessage("Boundaries should not intersect."));
          id = boundary_id;
        }
    return id;
  }

} // namespace FractureHealing
