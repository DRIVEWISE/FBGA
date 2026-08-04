/*
(***********************************************************************************)
(*                                                                                 *)
(* The FBGA project                                                               *)
(*                                                                                 *)
(* Copyright (c) 2025, Mattia Piazza                                               *)
(*                                                                                 *)
(*    Mattia Piazza                                                                *)
(*    Department of Industrial Engineering                                         *)
(*    University of Trento                                                         *)
(*    e-mail: mattia.piazza@unitn.it                                               *)
(*                                                                                 *)
(***********************************************************************************)
*/

#pragma once

#include <utils/types.hh>
#include <vector>
#include <string>
#include <cmath>
#include <limits>
#define STD_TOL 1.0e-10
#define STD_MAX_ITER 200
#define STD_VERBOSE "zero"

namespace fb::solvers
{

  using fb::utils::real;
  using fb::utils::integer;
  using fb::utils::QUIET_NAN;

  class BrentDekker
  {
  private:
    real tol = STD_TOL;
    integer max_iter = STD_MAX_ITER;
    std::string verbose = "iter";
    integer iter = 0;
    bool flag = false;
    real x_last = QUIET_NAN;
    std::vector<real> x_iter;

  public:
    BrentDekker() = default;
    explicit BrentDekker(real tol);
    BrentDekker(real tol, integer max_iter);
    BrentDekker(real tol, integer max_iter, const std::string &verbose);

    void set_tolerance(const real new_tol){ tol = new_tol; }
    void set_max_iter(const integer new_max_iter) { max_iter = new_max_iter; }
    void set_verbose(const std::string &new_verbose) { verbose = new_verbose; }
    void setup(real tol, integer max_iter, const std::string &verbose);

    // Templated on the callable so lambdas fully inline at the call site instead of paying
    // std::function's type-erasure overhead -- this is called per Brent-Dekker iteration on
    // every root-find in the FBGA/FBGA3D hot path. See FBGA3D_INTEGRATION_PLAN.md.
    template <class F>
    bool solve(F const &f, real a, real b, real &x0)
    {
      this->flag = false;
      real x1 = a;
      real f1 = f(x1);
      if (f1 == 0)
      {
        x0 = x1;
        this->flag = true;
        return true;
      }
      real x2 = b;
      real f2 = f(x2);
      if (f2 == 0)
      {
        x0 = x2;
        this->flag = true;
        return true;
      }
      if (f1 * f2 > 0.0)
      {
        x0 = std::numeric_limits<real>::quiet_NaN();
        return false;
      }
      real x3 = (x1 + x2) / static_cast<real>(2);
      for (integer ith = 0; ith <= this->max_iter; ith++)
      {
        const real f3 = f(x3);
        if (std::abs(f3) < this->tol)
        {
          x0 = x3;
          this->flag = true;
          break;
        }
        if (f1 * f3 < 0.0)
        {
          b = x3;
        }
        else
        {
          a = x3;
        }
        if (b - a < this->tol * std::max(std::abs(b), 1.0))
        {
          x0 = (x1 + x2) / static_cast<real>(2);
          this->flag = true;
          break;
        }
        const real denom = (f2 - f1) * (f3 - f1) * (f2 - f3);
        const real numer = x3 * (f1 - f2) * (f2 - f3 + f1) + f2 * x1 * (f2 - f3) + f1 * x2 * (f3 - f1);
        real dx = denom != 0.0 ? f3 * numer / denom : b - a;
        real x = x3 + dx;
        if ((b - x) * (x - a) < 0.0)
        {
          dx = (b - a) / static_cast<real>(2);
          x = a + dx;
        }
        if (x1 < x3)
        {
          x2 = x3;
          f2 = f3;
        }
        else
        {
          x1 = x3;
          f1 = f3;
        }
        if (std::abs(x - x3) < this->tol)
        {
          x0 = x;
          this->flag = true;
          break;
        }
        x3 = x;
      }
      this->iter = this->max_iter;
      this->x_last = x0;
      return this->flag;
    }
 };

}

