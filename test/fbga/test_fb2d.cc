#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <cmath>

#include <fbga/fb2d.hh>

using namespace fb::utils;
using namespace fb::fbga;
using Catch::Approx;

namespace
{
  // A symmetric friction-circle envelope (as in examples/fbga/basic.cc), with
  // longitudinal bound g*mu_x independent of ay when ay == 0.
  struct FrictionCircle
  {
    real mu_x = 1.0;
    real mu_y = 1.0;
    real g = 10.0;

    [[nodiscard]] real max_ax() const { return g * mu_x; }

    [[nodiscard]] real upper(real ay, real /*v*/) const
    {
      const real ay_norm = ay / g;
      return g * std::sqrt(-ay_norm * ay_norm + mu_y * mu_y) * mu_x / mu_y;
    }
    [[nodiscard]] real lower(real ay, real v) const { return -upper(ay, v); }
  };

  Fb2d make_straight_track_solver(const FrictionCircle &fc)
  {
    auto upper = [fc](real ay, real v) -> real { return fc.upper(ay, v); };
    auto lower = [fc](real ay, real v) -> real { return fc.lower(ay, v); };
    auto range_min = [fc](real /*v*/) -> real { return -fc.mu_y * fc.g; };
    auto range_max = [fc](real /*v*/) -> real { return fc.mu_y * fc.g; };

    return Fb2d(upper, lower, GgRangeMaxMin{range_min, range_max});
  }
}

TEST_CASE("Fb2d accelerates at the constant friction-circle limit on a straight track", "[fbga][Fb2d]")
{
  const FrictionCircle fc;
  Fb2d solver = make_straight_track_solver(fc);

  const std::vector<real> SS = {0.0, 10.0, 20.0, 30.0, 40.0};
  const std::vector<real> KK = {0.0, 0.0, 0.0, 0.0, 0.0}; // straight: zero curvature everywhere
  const real v0 = 20.0;
  const real total_length = SS.back() - SS.front();

  const real total_time = solver.compute(SS, KK, v0);

  // On a straight track the friction circle degenerates to a constant longitudinal
  // bound (mu_x * g), so the whole profile is a textbook constant-acceleration ramp.
  const real a = fc.max_ax();
  const real expected_vf = std::sqrt(v0 * v0 + 2.0 * a * total_length);
  const real expected_time = (expected_vf - v0) / a;

  CHECK(solver.get_dump().empty()); // no segment should have failed to solve

  CHECK(solver.evalV(0.0) == Approx(v0).margin(1e-6));
  CHECK(solver.evalV(total_length) == Approx(expected_vf).margin(1e-6));
  CHECK(total_time == Approx(expected_time).margin(1e-6));

  // Zero curvature -> zero lateral acceleration everywhere.
  CHECK(solver.evalAy(10.0) == Approx(0.0).margin(1e-9));
  // Constant longitudinal acceleration at the friction limit throughout.
  CHECK(solver.evalAx(10.0) == Approx(a).margin(1e-6));
  CHECK(solver.evalAx(30.0) == Approx(a).margin(1e-6));
}

TEST_CASE("Fb2d::is_in_range reports whether a point is within the g-g envelope", "[fbga][Fb2d]")
{
  const FrictionCircle fc;
  Fb2d solver = make_straight_track_solver(fc);

  // Inside the friction circle at v = 20: ax = 0, ay = 0 is trivially within bounds.
  CHECK(solver.is_in_range(0.0, 0.0, 20.0));
  // Far outside any physical friction circle.
  CHECK_FALSE(solver.is_in_range(1000.0, 0.0, 20.0));
}
