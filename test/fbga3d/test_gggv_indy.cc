#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <fbga3d/gggv_indy.hh>

using namespace fb::fbga3d;
using Catch::Approx;

// GggvIndy::setup_std() (invoked by the default constructor) loads spline data from
// "./data/INDY/*.npy" relative to the process working directory -- see this target's
// WORKING_DIRECTORY in test/fbga3d/CMakeLists.txt.

TEST_CASE("GggvIndy loads its default spline data without throwing", "[fbga3d][gggv_indy]")
{
  REQUIRE_NOTHROW(GggvIndy{});
}

TEST_CASE("GggvIndy::a_x_max/a_x_min bracket zero acceleration in opposite directions", "[fbga3d][gggv_indy]")
{
  GggvIndy gggv;
  const real V = 30.0;
  const real az_tilde = 10.0;
  CHECK(gggv.a_x_max(V, az_tilde) > 0.0);
  CHECK(gggv.a_x_min(V, az_tilde) < 0.0);
}

TEST_CASE("GggvIndy::a_x_push saturates to a_x_eng once engine power is the binding constraint", "[fbga3d][gggv_indy]")
{
  GggvIndy gggv;
  const real V = 30.0;
  const real az_tilde = 10.0;
  CHECK(gggv.a_x_push(0.0, V, az_tilde) <= gggv.a_x_eng(V) + 1e-9);
}

TEST_CASE("GggvIndy::set_scaling_factors rescales a_x_max linearly", "[fbga3d][gggv_indy]")
{
  GggvIndy gggv;
  const real V = 30.0;
  const real az_tilde = 10.0;
  const real baseline = gggv.a_x_max(V, az_tilde);

  ScalingGggvFactors scaling;
  scaling.ax_max_scale = 0.5;
  gggv.set_scaling_factors(scaling);

  CHECK(gggv.a_x_max(V, az_tilde) == Approx(baseline * 0.5));
}
