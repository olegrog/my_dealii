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

#include <deal.II/lac/vector.h>
#include <deal.II/numerics/vector_tools.h>

#include "Parameters.h"

// See details in Bailon-Plaza & Meulen J. theor. Biol. (2001) 212, 191-209

namespace
{
  template <int n_components>
  std::array<double, n_components>
  vector2array(const dealii::Vector<double> vec)
  {
    std::array<double, n_components> a;
    std::copy_n(vec.begin(), n_components, a.begin());
    return a;
  }
} // namespace

namespace FractureHealing
{
  Model::Model()
    : ParameterAcceptor("Model")
  {
    add_parameter("Haptotactic cell migration", Dh = 0.014);
    add_parameter("Haptotactic cell migration K", Kh = 0.25);
    add_parameter("Haptokinetic cell migration", Ck = 0.0034);
    add_parameter("Haptokinetic cell migration K", Kk = 0.5);
    add_parameter("Mesenchymal proliferation", Am0 = 1.01);
    add_parameter("Mesenchymal proliferation K", Km = 0.1);
    add_parameter("Chondrocyte proliferation", Ac0 = 1.01);
    add_parameter("Chondrocyte proliferation K", Kc = 0.1);
    add_parameter("Osteoblast proliferation", Ab0 = 0.202);
    add_parameter("Osteoblast proliferation K", Kb = 0.1);
    add_parameter("Mesenchymal limiting density", alpha_m = 1);
    add_parameter("Chondrocyte limiting density", alpha_c = 1);
    add_parameter("Osteoblast limiting density", alpha_b = 1);
    add_parameter("Diffentiation into chondrocytes", Y1 = 10);
    add_parameter("Diffentiation into chondrocytes H", H1 = 0.1);
    add_parameter("Diffentiation into osteoblasts", Y2 = 50);
    add_parameter("Diffentiation into osteoblasts H", H2 = 0.2);
    add_parameter("Endochondral ossification", Y3 = 1000);
    add_parameter("Endochondral ossification H", H3 = 0.1);
    add_parameter("Endochondral replacement B", Bec = 1.5);
    add_parameter("Bone decay", db = 0.1);
    add_parameter("Cartilage synthesis", Pcs = 0.2);
    add_parameter("Bone synthesis", Pbs = 2);
    add_parameter("Cartilage degradation 1", Qcd1 = 0.2);
    add_parameter("Cartilage degradation 2", Qcd2 = 2);
    add_parameter("Bone degradation", Qbd = 2);
    add_parameter("Chondrogenic GF diffusion", Dgc = 0.005);
    add_parameter("Osteogenic GF diffusion", Dgb = 0.005);
    add_parameter("Chondrogenic GF growth", Ggc = 50);
    add_parameter("Chondrogenic GF growth H", Hgc = 1);
    add_parameter("Chondrogenic GF growth K", Kgc = 0.1);
    add_parameter("Osteogenic GF growth", Ggb = 500);
    add_parameter("Osteogenic GF growth H", Hgb = 1);
    add_parameter("Chondrogenic GF decay", dgc = 100);
    add_parameter("Osteogenic GF decay", dgb = 100);
  }

  const std::vector<std::string> Model::component_names = {
    {"Mesenchymal cellular density",
     "Chondrocytes cellular density",
     "Osteoblasts cellular density",
     "Cartilage ECM density",
     "Bone ECM density",
     "Chondrogenic GF density",
     "Osteogenic GF density"}};



  void
  Model::get_diffusion_matrices(std::vector<matrix_type> &         matrices,
                                const std::vector<Vector<double>> &values) const
  {
    using std::pow;

    AssertDimension(matrices.size(), values.size());

    for (unsigned int i = 0; i < values.size(); i++)
      {
        matrix_type &A = matrices[i] = 0;

        const auto [cm, cc, cb, mc, mb, gc, gb] =
          vector2array<n_components>(values[i]);

        const double m = mc + mb;
        const double D = Dh * m / (pow(Kh, 2) + pow(m, 2));
        const double C = Ck / pow(Kk + m, 2);

        A[0][0] = D;
        A[0][3] = C * cm;
        A[5][5] = Dgc;
        A[6][6] = Dgb;
      }
  }



  void
  Model::get_reaction_matrices(std::vector<matrix_type> &         matrices,
                               const std::vector<Vector<double>> &values) const
  {
    using std::pow;

    AssertDimension(matrices.size(), values.size());

    for (unsigned int i = 0; i < values.size(); i++)
      {
        matrix_type &A = matrices[i] = 0;

        const auto [cm, cc, cb, mc, mb, gc, gb] =
          vector2array<n_components>(values[i]);

        const double m  = mc + mb;
        const double Am = Am0 * m / (pow(Km, 2) + pow(m, 2));
        const double Ac = Ac0 * m / (pow(Kc, 2) + pow(m, 2));
        const double Ab = Ab0 * m / (pow(Kb, 2) + pow(m, 2));
        const double F1 = Y1 * gb / (H1 + gb);
        const double F2 = Y2 * gc / (H2 + gc);
        const double F3 =
          pow(mc, 6) / (pow(Bec, 6) + pow(mc, 6)) * Y3 * gb / (H3 + gb);
        const double Egc =
          Ggc * gc / (Hgc + gc) * m / (pow(Kgc, 3) + pow(m, 3));
        const double Egb = Ggb * gb / (Hgb + gb);

        A[0][0] = Am * (1 - alpha_m * cm) - F1 - F2;
        A[1][1] = Ac * (1 - alpha_c * cc) - F3;
        A[1][0] = F2;
        A[2][2] = Ab * (1 - alpha_b * cb) - db;
        A[2][0] = F1;
        A[2][1] = F3;
        A[3][0] = (Pcs - Qcd1 * mc);
        A[3][1] = A[3][0];
        A[3][2] = -Qcd2 * mc;
        A[4][2] = (Pbs - Qbd * mb);
        A[5][5] = -dgc;
        A[5][1] = Egc;
        A[6][6] = -dgb;
        A[6][2] = Egb;
      }
  }



  void
  Model::get_source_vectors(std::vector<vector_type> &         vectors,
                            const std::vector<Vector<double>> &values) const
  {
    AssertDimension(vectors.size(), values.size());

    for (unsigned int i = 0; i < values.size(); i++)
      vectors[i] = 0;
  }



  template <int dim>
  Geometry<dim>::Geometry()
    : ParameterAcceptor("Geometry")
  {
    add_parameter("Size", size);
    add_parameter("Repetitions", repetitions);
    add_parameter("Cells to remove", n_cells_to_remove);
  }



  template <int dim>
  InitialValues<dim>::InitialValues()
    : ParameterAcceptor("Initial values")
    , FunctionParser<dim>(Model::n_components)
  {
    add_parameter("Project functions", project_functions = false);
  }



  template <int dim>
  void InitialValues<dim>::declare_parameters(ParameterHandler &prm)
  {
    for (auto component : Model::component_names)
      prm.declare_entry(component, "0.0");
  }



  template <int dim>
  void InitialValues<dim>::parse_parameters(ParameterHandler &prm)
  {
    std::vector<std::string> expressions;

    for (auto component : Model::component_names)
      expressions.push_back(prm.get(component));

    FunctionParser<dim>::initialize(
      FunctionParser<dim>::default_variable_names(),
      expressions,
      typename FunctionParser<dim>::ConstMap());
  }



  template <int dim>
  BoundaryValues<dim>::BoundaryValues()
    : ParameterAcceptor("Boundary values")
  {}



  template <int dim>
  void BoundaryValues<dim>::declare_parameters(ParameterHandler &prm)
  {
    for (unsigned int b = 0; b < max_n_boundaries; ++b)
      {
        prm.enter_subsection(subsection(b));
        {
          prm.declare_entry("Location", "0");
          for (auto component : Model::component_names)
            prm.declare_entry(component, "none");
        }
        prm.leave_subsection();
      }
  }



  template <int dim>
  void BoundaryValues<dim>::parse_parameters(ParameterHandler &prm)
  {
    for (unsigned int b = 0; b < max_n_boundaries; ++b)
      {
        prm.enter_subsection(subsection(b));

        std::vector<std::string> expressions;
        std::vector<bool>        cmask;

        for (auto component : Model::component_names)
          {
            const std::string expression = prm.get(component);
            if (expression == "none")
              {
                // Zero gradient BC
                expressions.push_back("-1");
                cmask.push_back(false);
              }
            else
              {
                // Dirichlet BC
                expressions.push_back(expression);
                cmask.push_back(true);
              }
          }

        if (std::find(begin(cmask), end(cmask), true) != end(cmask))
          boundaries.emplace(std::piecewise_construct,
                             std::forward_as_tuple(b),
                             std::forward_as_tuple(prm.get("Location"),
                                                   expressions,
                                                   cmask));
        prm.leave_subsection();
      }
  }



  template <int dim>
  void BoundaryValues<dim>::interpolate_boundary_values(
    const DoFHandler<dim> &                    dof_handler,
    std::map<types::global_dof_index, double> &boundary_values) const
  {
    for (const auto &[boundary_id, boundary] : boundaries)
      VectorTools::interpolate_boundary_values(dof_handler,
                                               boundary_id,
                                               boundary.values,
                                               boundary_values,
                                               boundary.component_mask);
  }



  template <int dim>
  BoundaryValues<dim>::Boundary::Boundary(
    const std::string &             location_expression,
    const std::vector<std::string> &value_expressions,
    const std::vector<bool> &       cmask)
    : location(location_expression,
               "",
               FunctionParser<dim>::default_variable_names())
    , values(Model::n_components)
    , component_mask(cmask)
  {
    values.initialize(FunctionParser<dim>::default_variable_names(),
                      value_expressions,
                      typename FunctionParser<dim>::ConstMap());
  }



  FiniteElements::FiniteElements()
    : ParameterAcceptor("Finite elements")
  {
    add_parameter("Polynomial degree", poly_degree = 1);
    add_parameter("Renumbering", renumbering = false);
  }



  MeshRefinement::MeshRefinement()
    : ParameterAcceptor("Mesh refinement")
  {
    add_parameter("Minimum level", j_min = 3);
    add_parameter("Maximum level", j_max = 5);
    add_parameter("Number of time steps", n_steps = 1);
    add_parameter("IC refining threshold", ic_upper = 0.99);
    add_parameter("Refining threshold", upper = 0.5);
    add_parameter("Coarsening threshold", lower = 0.1);
    add_parameter("Maximum number of cells", max_n_cells = 5000);
  }



  LinearSolver::LinearSolver()
    : ParameterAcceptor("Linear solver")
  {
    add_parameter("Solver name",
                  solver_name = "GMRES",
                  "",
                  prm,
                  Patterns::Selection("CG|BiCGStab|GMRES"));
    add_parameter("Max iterations", max_iter = 100);
    add_parameter("Tolerance", tol = 1e-8);
    add_parameter("Reduce", reduce = 1e-2);
    add_parameter("Preconditioner relaxation", preconditioner_relax = 1.0);
  }



  NonlinearSolver::NonlinearSolver()
    : ParameterAcceptor("Nonlinear solver")
  {
    add_parameter("Max iterations", max_iter = 15);
    add_parameter("Tolerance", tol = 1e-7);
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



  template class Geometry<2>;
  template class InitialValues<2>;
  template class BoundaryValues<2>;

} // namespace FractureHealing
