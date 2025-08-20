// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <prismspf/config.h>

#include "prismspf/core/types.h"

PRISMS_PF_BEGIN_NAMESPACE

class SimulationTime
{
public:
  /**
   * @brief Reset the simulation time and increment to 0.
   */
  void
  reset()
  {
    _time      = 0.0;
    _increment = 0;
  }

  /**
   * @brief Get the current simulation time.
   */
  [[nodiscard]] double
  current_time() const
  {
    return curr_time;
  }

  /**
   * @brief Get the current increment.
   */
  [[nodiscard]] unsigned int
  current_increment() const
  {
    return curr_increment;
  }

  /**
   * @brief Increment the simulation time by the given time step.
   */
  void
  increment(double time_step)
  {
    curr_time += time_step;
    ++curr_increment;
  }

private:
  /**
   * @brief The current simulation time.
   */
  double curr_time = 0.0;

  /**
   * @brief The time increment for the next simulation step.
   */
  unsigned int curr_increment = 0;
};

PRISMS_PF_END_NAMESPACE