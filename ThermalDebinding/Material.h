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

  class Material : ParameterAcceptor
  {
  public:
    struct PolymerSpecie
    {
      PolymerSpecie(double initialValue,
                    double constantRate,
                    double activationEnergy)
        : y(initialValue)
        , k(constantRate)
        , E(activationEnergy)
      {}
      double       y;
      const double k, E;
    };

    Material();

    void parse_parameters(ParameterHandler &prm) override;

    //- Return the polymer concentration change rate [1/s]
    double dydt(double T) const;

    //- Return the monomer pressure [Pa]
    double P(double rho, double T) const;

    //- Return the diffusion coefficient [m^2/s]
    double D1(double T) const;

    //- Return the mass diffusivity due to convection [m^2/s]
    double D2(double P) const;

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

    //- Return the volume fraction of pores
    double poreVolumeFraction() const
    {
      return 1 - ceramicVolumeFraction_ - polymerVolumeFraction();
    }

    //- Return the species
    const std::vector<PolymerSpecie> &species() const
    {
      return species_;
    }

    //- Integrate polymer specie concentrations over time
    void evolve(double T, double delta_t);

  private:
    //- Return the polymer specie concentration change rate [1/s]
    double dydt(const PolymerSpecie &specie, double y, double T) const;

    //- Return the dynamic viscosity [kg/m/s]
    double mu() const
    {
      return mu_;
    }

    //- Return the permeability coefficient [m^2]
    double K() const;

    //- Return the volume fraction of polymers
    double polymerVolumeFraction() const;

    //- Density of the polymer [kg/m^3]
    double polymerRho_;

    //- Molar mass of the monomer [kg/mol]
    double monomerW_;

    //- Volume fraction of ceramic inclusions []
    double ceramicVolumeFraction_;

    //- Initial porosity []
    double initialPorosity_;

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

    //- Polymer species
    std::vector<PolymerSpecie> species_;

    //- Temporary array for initial values of polymer species
    std::vector<double> y0_;

    //- Temporary array for constant rates of polymer species
    std::vector<double> k_;

    //- Temporary array for activation energies of polymer species
    std::vector<double> E_;
  };

} // namespace ThermalDebinding
