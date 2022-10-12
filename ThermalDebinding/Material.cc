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

#include <numeric>

#include <deal.II/sundials/arkode.h>

#include "Material.h"

namespace ThermalDebinding
{
  Material::Material()
    : ParameterAcceptor("Material")
  {
    add_parameter("Polymer density", polymerRho_ = 1e3);
    add_parameter("Monomer molar mass", monomerW_ = 100e-3);
    add_parameter("Ceramic volume fraction", ceramicVolumeFraction_ = 0.63);
    add_parameter("Initial porosity", initialPorosity_ = 0.03);
    add_parameter("Dynamic viscosity", mu_ = 2e-3);
    add_parameter("Mean particle size", meanParticleSize_ = 1e-7);
    add_parameter("Particle size exponent", particleSizeExponent_ = 2);

    enter_subsection("Polymer composition");
    {
      add_parameter("Initial mass fraction", y0_);
      add_parameter("Constant rate", k_);
      add_parameter("Activation energy", E_);
    }
    leave_subsection();

    enter_subsection("Diffusion");
    {
      add_parameter("Pre-exponential factor", D0_ = 6.92e-4);
      add_parameter("Activation energy", Ea_ = 38.37e3);
    }
    leave_subsection();
  }


  // Public member functions


  void Material::parse_parameters(ParameterHandler &)
  {
    Assert(y0_.size() == k_.size(),
           ExcDimensionMismatch(y0_.size(), k_.size()));
    Assert(y0_.size() == E_.size(),
           ExcDimensionMismatch(y0_.size(), E_.size()));

    for (unsigned i = 0; i < y0_.size(); i++)
      species_.emplace_back(PolymerSpecie(y0_[i], k_[i], E_[i]));
  }



  double Material::dydt(double T) const
  {
    return std::accumulate(species_.cbegin(),
                           species_.cend(),
                           0.,
                           [this, T](double res, const PolymerSpecie &specie) {
                             return res + dydt(specie, specie.y, T);
                           });
  }



  double Material::P(double rho, double T) const
  {
    const double R     = 8.3145; // J/mol/K
    const double small = 1e-16;
    return rho * R * T / monomerW_ / (poreVolumeFraction() + small);
  }



  double Material::D1(double T) const
  {
    const double R = 8.3145; // J/mol/K
    return D0_ * exp(-Ea_ / R / T) * std::pow(1 - ceramicVolumeFraction_, 1.5);
  }



  double Material::D2(double P) const
  {
    return K() / mu() / poreVolumeFraction() * P;
  }



  void Material::evolve(double T, double delta_t)
  {
    using VectorType = Vector<double>;

    SUNDIALS::ARKode<VectorType>::AdditionalData params(
      /*initial_time = */ 0, /*final_time = */ delta_t);
    SUNDIALS::ARKode<VectorType> ode(params);

    ode.explicit_function =
      [this, T](double, const VectorType &y, VectorType &ydot) -> int {
      for (unsigned i = 0; i < y.size(); i++)
        ydot[i] = dydt(species_[i], y[i], T);
      return 0;
    };

    VectorType y(species_.size());
    for (unsigned i = 0; i < y.size(); i++)
      y[i] = species_[i].y;

    ode.solve_ode(y);

    for (unsigned i = 0; i < y.size(); i++)
      species_[i].y = y[i];
  }


  // Private member functions


  double Material::dydt(const PolymerSpecie &specie, double y, double T) const
  {
    const double R = 8.3145; // J/mol/K
    return -specie.k * y * std::exp(-specie.E / R / T);
  }



  double Material::K() const
  {
    const double K0   = std::pow(meanParticleSize_, 2) / 180;
    const double phiM = poreVolumeFraction();

    return K0 *
           std::pow(phiM / (1 - ceramicVolumeFraction_),
                    particleSizeExponent_) *
           std::pow(phiM, 3) / std::pow(1 - phiM, 2);
  }



  double Material::polymerVolumeFraction() const
  {
    double y = std::accumulate(species_.cbegin(),
                               species_.cend(),
                               0,
                               [](double res, const PolymerSpecie &specie) {
                                 return res + specie.y;
                               });
    return y * initialPolymerFraction();
  }

} // namespace ThermalDebinding
