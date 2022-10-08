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

#include <deal.II/base/parameter_handler.h>

namespace ThermalDebinding
{
  using namespace dealii;

  struct Problem
  {
    double size;
    double T0;
    double heating_rate;

    static void declare_parameters(ParameterHandler &prm);

    void parse_parameters(ParameterHandler &prm);
  };

  class Material
  {
  public:
    void        parse_parameters(ParameterHandler &prm);
    static void declare_parameters(ParameterHandler &prm);

    //- Return the polymer concentration w.r.t. the initial one
    double y(double time) const
    {
      return y0_ * exp(-degradationRate_ * time);
    }

    //- Return the polymer concentration change rate [1/s]
    double dydt(double time) const
    {
      return -degradationRate_ * y(time);
    }

    //- Return the polymer density [kg/m^3]
    double polymerRho() const
    {
      return polymerRho_;
    }

    //- Return the initial volume fraction of polymers
    double initialPolymerFraction() const
    {
      return 1 - ceramicVolumeFraction_ - initialPorosity_;
    }

    //- Return the volume fraction of polymers
    double polymerVolumeFraction(double time) const
    {
      return y(time) * initialPolymerFraction();
    }

    //- Return the volume fraction of pores
    double poreVolumeFraction(double time) const
    {
      return 1 - ceramicVolumeFraction_ - polymerVolumeFraction(time);
    }

    //- Return the monomer pressure [Pa]
    double P(double time, double rho, double T) const
    {
      const double R     = 8.3145; // J/mol/K
      const double small = 1e-16;
      return rho * R * T / monomerW_ / (poreVolumeFraction(time) + small);
    }

    //- Return the diffusion coefficient [m^2/s]
    double D(double T) const
    {
      const double R = 8.3145; // J/mol/K
      return D0_ * exp(-Ea_ / R / T) *
             std::pow(1 - ceramicVolumeFraction_, 1.5);
    }

    //- Return the dynamic viscosity [kg/m/s]
    double mu() const
    {
      return mu_;
    }

    //- Return the permeability coefficient [m^2]
    double K(double time) const
    {
      const double K0   = std::pow(meanParticleSize_, 2) / 180;
      const double phiM = poreVolumeFraction(time);

      return K0 *
             std::pow(phiM / (1 - ceramicVolumeFraction_),
                      particleSizeExponent_) *
             std::pow(phiM, 3) / std::pow(1 - phiM, 2);
    }

  private:
    //- Density of the polymer [kg/m^3]
    double polymerRho_;

    //- Molar mass of the monomer [kg/mol]
    double monomerW_;

    //- Volume fraction of ceramic inclusions []
    double ceramicVolumeFraction_;

    //- Initial porosity []
    double initialPorosity_;

    //- Initial polymer concentration []
    double y0_;

    //- Thermal degradation rate of the polymer [1/s]
    double degradationRate_;

    //- Bulk diffusion coefficient [m^2/s]
    double D0_;

    //- Diffusion activation energy [J/mol]
    double Ea_;

    //- Dynamic viscosity [kg/m/s]
    double mu_;

    //- Mean ceramic particle size [m]
    double meanParticleSize_;

    //- Model parameter for permeability []
    double particleSizeExponent_;
  };

  class Time
  {
  public:
    Time()
      : step_(0)
      , current_(0.0)
    {}

    static void declare_parameters(ParameterHandler &prm);

    void parse_parameters(ParameterHandler &prm);


    double operator()() const
    {
      return current_;
    }
    double end() const
    {
      return end_;
    }
    double delta() const
    {
      return delta_;
    }
    double theta() const
    {
      return theta_;
    }
    unsigned int step() const
    {
      return step_;
    }
    void operator++()
    {
      current_ += delta_;
      ++step_;
    }
    bool loop()
    {
      operator++();
      return current_ <= end_;
    }

  private:
    unsigned int step_;
    double       current_;
    double       end_;
    double       delta_;
    double       theta_;
  };

  struct FiniteElements
  {
    unsigned int poly_degree;
    unsigned int quad_order;

    static void declare_parameters(ParameterHandler &prm);

    void parse_parameters(ParameterHandler &prm);
  };

  struct MeshRefinement
  {
    unsigned int j_min;
    unsigned int j_max;
    unsigned int n_steps;
    double       upper;
    double       lower;
    double       max_n_cells;

    static void declare_parameters(ParameterHandler &prm);

    void parse_parameters(ParameterHandler &prm);
  };

  struct LinearSolver
  {
    unsigned int max_iter;
    double       tol;
    double       reduce;
    double       preconditioner_relax;

    static void declare_parameters(ParameterHandler &prm);

    void parse_parameters(ParameterHandler &prm);
  };

  struct NonlinearSolver
  {
    unsigned int max_iter;
    double       tol;

    static void declare_parameters(ParameterHandler &prm);

    void parse_parameters(ParameterHandler &prm);
  };

  struct Output
  {
    bool         write_vtk_files;
    unsigned int n_steps;
    unsigned int verbosity;

    static void declare_parameters(ParameterHandler &prm);

    void parse_parameters(ParameterHandler &prm);
  };

  struct Parameters
  {
    Problem         problem;
    Material        material;
    Time            time;
    FiniteElements  fe;
    MeshRefinement  mr;
    LinearSolver    ls;
    NonlinearSolver ns;
    Output          output;

    Parameters(const std::string &input_file);
    static void declare_parameters(ParameterHandler &prm);
    void        parse_parameters(ParameterHandler &prm);
  };

} // namespace ThermalDebinding
