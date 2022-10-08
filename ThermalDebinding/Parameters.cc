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
  Parameters::Parameters(const std::string &input_file)
  {
    ParameterHandler prm;
    declare_parameters(prm);
    prm.parse_input(input_file);
    parse_parameters(prm);
  }

  void Parameters::declare_parameters(ParameterHandler &prm)
  {
    Problem::declare_parameters(prm);
    Material::declare_parameters(prm);
    Time::declare_parameters(prm);
    FiniteElements::declare_parameters(prm);
    MeshRefinement::declare_parameters(prm);
    LinearSolver::declare_parameters(prm);
    NonlinearSolver::declare_parameters(prm);
    Output::declare_parameters(prm);
  }

  void Parameters::parse_parameters(ParameterHandler &prm)
  {
    problem.parse_parameters(prm);
    material.parse_parameters(prm);
    time.parse_parameters(prm);
    fe.parse_parameters(prm);
    mr.parse_parameters(prm);
    ls.parse_parameters(prm);
    ns.parse_parameters(prm);
    output.parse_parameters(prm);
  }



  void Problem::declare_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("Problem");
    {
      prm.declare_entry("Hyper cube size", "1e-2", Patterns::Double(0));
      prm.declare_entry("Initial temperature", "300", Patterns::Double(0));
      prm.declare_entry("Heating rate", "1.667e-3", Patterns::Double(0));
    }
    prm.leave_subsection();
  }

  void Problem::parse_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("Problem");
    {
      size         = prm.get_double("Hyper cube size");
      T0           = prm.get_double("Initial temperature");
      heating_rate = prm.get_double("Heating rate");
    }
    prm.leave_subsection();
  }



  void Material::declare_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("Material");
    {
      prm.declare_entry("Polymer density", "1e3", Patterns::Double(0));
      prm.declare_entry("Monomer molar mass", "1e-1", Patterns::Double(0));
      prm.declare_entry("Ceramic volume fraction",
                        "0.63",
                        Patterns::Double(0, 1));
      prm.declare_entry("Initial porosity", "0.03", Patterns::Double(0, 1));
      prm.declare_entry("Dynamic viscosity", "2e-3", Patterns::Double(0));
      prm.declare_entry("Mean particle size", "1e-7", Patterns::Double(0));
      prm.declare_entry("Particle size exponent", "2", Patterns::Double(0));
      prm.enter_subsection("Polymer composition");
      {
        prm.declare_entry("Number of components", "1", Patterns::Integer(0));
        prm.enter_subsection("Polymer1");
        {
          prm.declare_entry("Initial mass fraction",
                            "1",
                            Patterns::Double(0, 1));
          prm.declare_entry("Degradation rate", "1e-5", Patterns::Double(0));
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
      prm.enter_subsection("Diffusion");
      {
        prm.declare_entry("Pre-exponential factor",
                          "6.92e-4",
                          Patterns::Double(0));
        prm.declare_entry("Activation energy", "38.37e3", Patterns::Double(0));
      }
      prm.leave_subsection();
    }
    prm.leave_subsection();
  }

  void Material::parse_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("Material");
    {
      polymerRho_            = prm.get_double("Polymer density");
      monomerW_              = prm.get_double("Monomer molar mass");
      ceramicVolumeFraction_ = prm.get_double("Ceramic volume fraction");
      initialPorosity_       = prm.get_double("Initial porosity");
      mu_                    = prm.get_double("Dynamic viscosity");
      meanParticleSize_      = prm.get_double("Mean particle size");
      particleSizeExponent_  = prm.get_double("Particle size exponent");

      prm.enter_subsection("Polymer composition");
      {
        prm.enter_subsection("Polymer1");
        {
          y0_              = prm.get_double("Initial mass fraction");
          degradationRate_ = prm.get_double("Degradation rate");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();

      prm.enter_subsection("Diffusion");
      {
        D0_ = prm.get_double("Pre-exponential factor");
        Ea_ = prm.get_double("Activation energy");
      }
      prm.leave_subsection();
    }
    prm.leave_subsection();
  }



  void Time::declare_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("Time");
    {
      prm.declare_entry("End time", "500e3", Patterns::Double(1e-100));
      prm.declare_entry("Step size", "1e3", Patterns::Double(1e-100));
      prm.declare_entry("Theta", "1", Patterns::Double(0, 1));
    }
    prm.leave_subsection();
  }

  void Time::parse_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("Time");
    {
      end_   = prm.get_double("End time");
      delta_ = prm.get_double("Step size");
      theta_ = prm.get_double("Theta");
    }

    prm.leave_subsection();
  }



  void FiniteElements::declare_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("Finite elements");
    {
      prm.declare_entry("Polynomial degree", "2", Patterns::Integer(0));
      prm.declare_entry("Quadrature order", "3", Patterns::Integer(0));
    }
    prm.leave_subsection();
  }

  void FiniteElements::parse_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("Finite elements");
    {
      poly_degree = prm.get_integer("Polynomial degree");
      quad_order  = prm.get_integer("Quadrature order");
    }
    prm.leave_subsection();
  }



  void MeshRefinement::declare_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("Mesh refinement");
    {
      prm.declare_entry("Minimum level", "4", Patterns::Integer(1));
      prm.declare_entry("Maximum level", "8", Patterns::Integer(1));
      prm.declare_entry("Number of time steps", "3", Patterns::Integer(1));
      prm.declare_entry("Refining threshold", "1", Patterns::Double(0, 1));
      prm.declare_entry("Coarsening threshold", "0", Patterns::Double(0, 1));
      prm.declare_entry("Maximum number of cells",
                        "1000",
                        Patterns::Integer(1));
    }
    prm.leave_subsection();
  }

  void MeshRefinement::parse_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("Mesh refinement");
    {
      j_min       = prm.get_integer("Minimum level");
      j_max       = prm.get_integer("Maximum level");
      n_steps     = prm.get_integer("Number of time steps");
      upper       = prm.get_double("Refining threshold");
      lower       = prm.get_double("Coarsening threshold");
      max_n_cells = prm.get_double("Maximum number of cells");
    }
    prm.leave_subsection();
  }



  void LinearSolver::declare_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("Linear solver");
    {
      prm.declare_entry("Max iterations", "15", Patterns::Integer(0));
      prm.declare_entry("Tolerance", "1e-8", Patterns::Double(0, 1));
      prm.declare_entry("Reduce", "1e-2", Patterns::Double(0, 1));
      prm.declare_entry("Preconditioner relaxation",
                        "1",
                        Patterns::Double(0, 2));
    }
    prm.leave_subsection();
  }

  void LinearSolver::parse_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("Linear solver");
    {
      max_iter             = prm.get_integer("Max iterations");
      tol                  = prm.get_double("Tolerance");
      reduce               = prm.get_double("Reduce");
      preconditioner_relax = prm.get_double("Preconditioner relaxation");
    }
    prm.leave_subsection();
  }



  void NonlinearSolver::declare_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("Nonlinear solver");
    {
      prm.declare_entry("Max iterations", "15", Patterns::Integer(0));
      prm.declare_entry("Tolerance", "1e-8", Patterns::Double(0.0));
    }
    prm.leave_subsection();
  }

  void NonlinearSolver::parse_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("Nonlinear solver");
    {
      max_iter = prm.get_integer("Max iterations");
      tol      = prm.get_double("Tolerance");
    }
    prm.leave_subsection();
  }



  void Output::declare_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("Output");
    {
      prm.declare_entry("Write VTK files", "true", Patterns::Bool());
      prm.declare_entry("Number of time steps", "1", Patterns::Integer(0));
      prm.declare_entry("Log verbosity", "0", Patterns::Integer(0));
    }
    prm.leave_subsection();
  }

  void Output::parse_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("Output");
    {
      write_vtk_files = prm.get_bool("Write VTK files");
      n_steps         = prm.get_integer("Number of time steps");
      verbosity       = prm.get_integer("Log verbosity");
    }
    prm.leave_subsection();
  }


} // namespace ThermalDebinding
