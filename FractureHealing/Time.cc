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
    , output_index_(0)
    , current_(0.0)
    , output_time_(false)
  {
    add_parameter("End time", end_ = 1);
    add_parameter("Step size", delta_ = 0.01);
    add_parameter("Theta", theta_ = 0.5);
    add_parameter("Time interval for output", time_for_output_ = 0.1);
    add_parameter("Iteration interval for output", n_steps_for_output_ = 1);
    add_parameter("Adaptive", adaptive_ = false);
    add_parameter("Tolerance", tol_ = 1e-2);
  }



  void Time::operator++()
  {
    current_ += delta_;
    step_++;

    output_time_ = false;
    if (time_for_output_ > 0)
      {
        const unsigned int index((current_ + delta_ / 2) / time_for_output_);
        if (index > output_index_)
          {
            output_index_ = index;
            output_time_  = true;
          }
      }
    else
      {
        output_index_ = step_ / n_steps_for_output_;
        output_time_  = step_ % n_steps_for_output_ == 0;
      }
  }



  void Time::update_delta(double residual_norm)
  {
    delta_ *= tol_ / residual_norm;

    // Adjust time step size to the output time interval
    if (time_for_output_ > 0)
      {
        const double time_till_output =
          std::max(0.0, (output_index_ + 1) * time_for_output_ - current_);
        const unsigned long n_steps =
          std::max(1L, std::lround(time_till_output / delta_));
        const double new_delta = time_till_output / n_steps;

        if (new_delta >= delta_)
          {
            delta_ = std::min(new_delta, 2.0 * delta_);
          }
        else
          {
            delta_ = std::max(new_delta, 0.2 * delta_);
          }
      }
  }

} // namespace FractureHealing
