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

} // namespace ThermalDebinding
