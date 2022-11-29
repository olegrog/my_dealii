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

namespace FractureHealing
{
  Time::Time()
    : ParameterAcceptor("Time")
    , step_(0)
    , current_(0.0)
  {
    add_parameter("End time", end_ = 1);
    add_parameter("Step size", delta_ = 0.1);
    add_parameter("Theta", theta_ = 0.5);
  }

} // namespace FractureHealing
