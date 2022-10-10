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

  class Material : ParameterAcceptor
  {
  public:
    Material();

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
    double D1(double T) const
    {
      const double R = 8.3145; // J/mol/K
      return D0_ * exp(-Ea_ / R / T) *
             std::pow(1 - ceramicVolumeFraction_, 1.5);
    }

    //- Return the mass diffusivity due to convection [m^2/s]
    double D2(double time, double P) const
    {
      return K(time) / mu() / poreVolumeFraction(time) * P;
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

  class Time : ParameterAcceptor
  {
  public:
    Time();

    //- Return the current time
    double operator()() const
    {
      return current_;
    }

    //- Return the end time
    double end() const
    {
      return end_;
    }

    //- Return the time step size
    double delta() const
    {
      return delta_;
    }

    //- Return parameter of the time-marching scheme
    double theta() const
    {
      return theta_;
    }

    //- Return the current iteration number
    unsigned int step() const
    {
      return step_;
    }

    //- Increment the time
    void operator++()
    {
      current_ += delta_;
      ++step_;
    }

    //- Increment the time and check if the end time is reached
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

  struct FiniteElements : ParameterAcceptor
  {
    FiniteElements();

    unsigned int poly_degree;
    unsigned int quad_order;
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

    bool         write_vtk_files;
    unsigned int n_steps;
    unsigned int verbosity;
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
  };

} // namespace ThermalDebinding
