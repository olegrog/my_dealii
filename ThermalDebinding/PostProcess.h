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
  using namespace dealii;

  template <int dim>
  class ComputePressure : public DataPostprocessorScalar<dim>
  {
  public:
    ComputePressure(const Material &material, double T)
      : DataPostprocessorScalar<dim>("pressure", update_values)
      , material_(material)
      , T_(T)
    {}

    virtual void evaluate_scalar_field(
      const DataPostprocessorInputs::Scalar<dim> &inputs,
      std::vector<Vector<double>> &computed_quantities) const override
    {
      for (unsigned int i = 0; i < computed_quantities.size(); ++i)
        computed_quantities[i] = material_.P(inputs.solution_values[i], T_);
    }

  private:
    const Material &material_;
    const double    T_;
  };

} // namespace ThermalDebinding
