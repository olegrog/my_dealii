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

namespace ThermalDebinding
{

  class Polymer
  {
  public:
    //- Constructor
    Polymer()
      : polymerRho_(1e3)
      , monomerW_(100e-3)
      , totalVolumeFraction_(0.37)
      , initialPorosity_(0.03)
      , y0_(1)
      , degradationRate_(1e-5)
      , D0_(6.92e-4)
      , Ea_(38.37e3)
    {}

    //- Return the polymer concentration
    double y(double time) const
    {
      return y0_*exp(-degradationRate_*time);
    }

    //- Return the polymer concentration change rate
    double dydt(double time) const
    {
      return -degradationRate_*y(time);
    }

    //- Return the polymer density [kg/m^3]
    double rho() const
    {
      return polymerRho_;
    }

    //- Return the initial volume fraction of polymer
    double initialVolumeFraction() const
    {
      return totalVolumeFraction_ - initialPorosity_;
    }

    //- Return the volume fraction of polymer
    double volumeFraction(double time) const
    {
      return y(time)*initialVolumeFraction();
    }

    //- Return the volume fraction of pores
    double poresFraction(double time) const
    {
      return totalVolumeFraction_ - volumeFraction(time);
    }

    //- Return the total volume fraction
    double totalVolumeFraction() const
    {
      return totalVolumeFraction_;
    }

    //- Return the monomer pressure [Pa]
    double pressure(double time, double rho, double T) const
    {
      const double R = 8.3145; // J/mol/K
      const double small = 1e-16;
      return rho*R*T/monomerW_/(poresFraction(time) + small);
    }

    //- Return the diffusion coefficient [m^2/s]
    double diffusion(double T) const
    {
      const double R = 8.3145; // J/mol/K
      return D0_*exp(-Ea_/R/T)*std::pow(totalVolumeFraction_, 1.5);
    }

  private:
        //- Density of the polymer [kg/m^3]
        const double polymerRho_;

        //- Molar mass of the monomer [kg/mol]
        const double monomerW_;

        //- Total volume fraction of polymer and pores []
        const double totalVolumeFraction_;

        //- Initial porosity []
        const double initialPorosity_;

        //- Initial polymer concentration []
        const double y0_;

        //- Thermal degradation rate of the polymer [1/s]
        const double degradationRate_;

        //- Bulk diffusion coefficient [m^2/s]
        const double D0_;

        //- Diffusion activation energy [J/mol]
        const double Ea_;
  };

} // namespace ThermalDebinding
