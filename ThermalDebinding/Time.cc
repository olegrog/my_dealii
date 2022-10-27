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

#include "Time.h"

namespace ThermalDebinding
{
  Time::Time()
    : ParameterAcceptor("Time")
    , step_(0)
    , current_(0.0)
  {
    add_parameter("End time", end_ = 500e3);
    add_parameter("Step size", delta_ = 1e3);
    add_parameter("Theta", theta_ = 0.5);
    add_parameter("Adaptive", adaptive_ = false);
    add_parameter("Maximum step size", max_delta_ = 10e3);
    add_parameter("Monomer production per step", delta_y_per_step_ = 1e-2);
  }

} // namespace ThermalDebinding
