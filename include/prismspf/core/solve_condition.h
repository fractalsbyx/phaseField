// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/base/exceptions.h>
#include <deal.II/matrix_free/evaluation_flags.h>

#include <prismspf/core/type_enums.h>
#include <prismspf/core/types.h>

#include <prismspf/config.h>

#include <functional>

PRISMS_PF_BEGIN_NAMESPACE

namespace SolveConditions
{
  struct SolveCondition
  {
    virtual bool
    should_solve(unsigned int increment);
    virtual ~SolveCondition();
  };

  struct Always : public SolveCondition
  {
    bool
    should_solve([[maybe_unused]] unsigned int increment) override
    {
      return true;
    }
  };

  struct Postprocess : public SolveCondition
  {
    bool
    should_solve([[maybe_unused]] unsigned int increment) override
    {
      return false;
    }
  };

  struct EveryN : public SolveCondition
  {
    explicit EveryN(unsigned int period)
      : n(period)
    {}

    bool
    should_solve(unsigned int increment) override
    {
      return !bool(increment % n);
    }

  private:
    unsigned int n;
  };

  // TODO: Add conditions for nucleation and AMR. This requires user_inputs to be done
  // nicely,
} // namespace SolveConditions

PRISMS_PF_END_NAMESPACE
