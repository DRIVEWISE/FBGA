#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <fbga3d/gggv_moto.hh>

using namespace fb::fbga3d;
using Catch::Approx;

TEST_CASE("GggvMoto::a_x_eng decreases with speed (power-limited)", "[fbga3d][gggv_moto]")
{
  GggvMoto gggv;
  const real a_low = gggv.a_x_eng(10.0);
  const real a_high = gggv.a_x_eng(50.0);
  CHECK(a_high < a_low);
}

TEST_CASE("GggvMoto::a_y_lim is symmetric adherence-only bound", "[fbga3d][gggv_moto]")
{
  GggvMoto gggv;
  const real az_tilde = 9.81;
  CHECK(gggv.a_y_lim(20.0, az_tilde) == Approx(az_tilde * MotoData{}.mu_Y));
}

TEST_CASE("GggvMoto::a_x_push/a_x_pull bracket zero lateral acceleration", "[fbga3d][gggv_moto]")
{
  GggvMoto gggv;
  const real V = 10.0;
  const real az_tilde = 9.81;
  const real a_x_max = gggv.a_x_push(0.0, V, az_tilde);
  const real a_x_min = gggv.a_x_pull(0.0, V, az_tilde);
  CHECK(a_x_max > 0.0);
  CHECK(a_x_min < 0.0);
  CHECK(a_x_max > a_x_min);
}

TEST_CASE("GggvMoto::a_x_push shrinks towards zero as lateral acceleration approaches its limit", "[fbga3d][gggv_moto]")
{
  GggvMoto gggv;
  const real V = 10.0;
  const real az_tilde = 9.81;
  const real a_y_lim = gggv.a_y_lim(V, az_tilde);
  const real a_x_center = gggv.a_x_push(0.0, V, az_tilde);
  const real a_x_at_limit = gggv.a_x_push(a_y_lim, V, az_tilde);
  CHECK(a_x_at_limit < a_x_center);
}
