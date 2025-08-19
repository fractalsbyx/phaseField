// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <prismspf/config.h>

#include "prismspf/core/types.h"

PRISMS_PF_BEGIN_NAMESPACE

struct SimulationTime
{
  /**
   * @brief The current simulation time.
   */
  double current_time = 0.0;

  /**
   * @brief The time increment for the next simulation step.
   */
  unsigned int current_increment = 0;

  /**
   * @brief Reset the simulation time and increment to 0.
   */
  void
  reset()
  {
    current_time      = 0.0;
    current_increment = 0;
  }

  /**
   * @brief Increment the simulation time by the given time step.
   */
  void
  increment(double time_step)
  {
    current_time += time_step;
    ++current_increment;
  }
};

PRISMS_PF_END_NAMESPACE