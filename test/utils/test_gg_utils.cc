#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <utils/gg_utils.hh>

using namespace fb::utils;
using Catch::Approx;

TEST_CASE("clip clamps a value into [a, b]", "[utils][clip]")
{
  CHECK(clip(0.5, -1.0, 1.0) == Approx(0.5));
  CHECK(clip(-5.0, -1.0, 1.0) == Approx(-1.0));
  CHECK(clip(5.0, -1.0, 1.0) == Approx(1.0));
}

TEST_CASE("piramid returns the max distance to the unit square boundary", "[utils][piramid]")
{
  CHECK(piramid(0.0, 0.0) == Approx(-1.0));
  CHECK(piramid(1.0, 0.0) == Approx(0.0));
  CHECK(piramid(2.0, 0.0) == Approx(1.0));
}

TEST_CASE("square_conversion normalizes a into [-1, 1] (clamped)", "[utils][square_conversion]")
{
  CHECK(square_conversion(0.0, -1.0, 1.0) == Approx(0.0));
  CHECK(square_conversion(1.0, -1.0, 1.0) == Approx(1.0));
  CHECK(square_conversion(-1.0, -1.0, 1.0) == Approx(-1.0));
}

TEST_CASE("signed_distance combines square_conversion and piramid", "[utils][signed_distance]")
{
  // center of the box on both axes -> deep inside -> negative distance
  CHECK(signed_distance(0.0, -1.0, 1.0, 0.0, -1.0, 1.0) < 0.0);
  // exactly on the ax boundary -> zero distance
  CHECK(signed_distance(1.0, -1.0, 1.0, 0.0, -1.0, 1.0) == Approx(0.0));
}

TEST_CASE("find_interval locates the bracketing segment", "[utils][find_interval]")
{
  const std::vector<real> vec = {0.0, 1.0, 2.0, 3.0, 4.0};

  CHECK(find_interval(vec, 0.5) == 0);
  CHECK(find_interval(vec, 2.5) == 2);
  CHECK(find_interval(vec, -1.0) == 0);       // below range -> clamps to first segment
  CHECK(find_interval(vec, 10.0) == vec.size() - 2); // above range -> clamps to last segment
}

TEST_CASE("find_interval_binary_search agrees with find_interval", "[utils][find_interval]")
{
  const std::vector<real> vec = {0.0, 1.0, 2.0, 3.0, 4.0};
  for (real v : {0.5, 1.5, 2.9})
  {
    CHECK(find_interval_binary_search(vec, v) == find_interval(vec, v));
  }
}

TEST_CASE("LinearInterpolator reproduces the sampled points and interpolates linearly", "[utils][LinearInterpolator]")
{
  const std::vector<real> X = {0.0, 1.0, 2.0};
  const std::vector<real> Y = {0.0, 2.0, 4.0}; // Y = 2*X

  LinearInterpolator interp(X, Y);

  CHECK(interp.eval(0.0) == Approx(0.0));
  CHECK(interp.eval(1.0) == Approx(2.0));
  CHECK(interp.eval(0.5) == Approx(1.0));
}

TEST_CASE("LinearInterpolatorSet interpolates the named column", "[utils][LinearInterpolatorSet]")
{
  const std::vector<real> X = {0.0, 1.0, 2.0};
  LinearInterpolatorSet set(X, { {0.0, 2.0, 4.0}, {0.0, -1.0, -2.0} }, {"double", "negative_half"});

  CHECK(set.eval("double", 0.5) == Approx(1.0));
  CHECK(set.eval("negative_half", 1.0) == Approx(-1.0));
}

TEST_CASE("BilinearInterpolator reproduces a bilinear grid", "[utils][BilinearInterpolator]")
{
  const std::vector<real> X = {0.0, 1.0};
  const std::vector<real> Y = {0.0, 1.0};
  // Z laid out row-major over X (see BilinearInterpolator::eval indexing): Z[i * ny + j]
  const std::vector<real> Z = {0.0, 1.0, 1.0, 2.0}; // Z(x,y) = x + y

  BilinearInterpolator interp(X, Y, Z);

  CHECK(interp.eval(0.0, 0.0) == Approx(0.0));
  CHECK(interp.eval(1.0, 1.0) == Approx(2.0));
  CHECK(interp.eval(0.5, 0.5) == Approx(1.0));
}

TEST_CASE("linspace generates N evenly spaced points", "[utils][linspace]")
{
  std::vector<real> vec;
  linspace(vec, 0.0, 10.0, 5);

  REQUIRE(vec.size() == 5);
  CHECK(vec.front() == Approx(0.0));
  CHECK(vec.back() == Approx(10.0));
  CHECK(vec[1] == Approx(2.5));
}

TEST_CASE("compute_finite_difference approximates the derivative of a linear function", "[utils][compute_finite_difference]")
{
  const std::vector<double> X = {0.0, 1.0, 2.0, 3.0};
  const std::vector<double> Y = {0.0, 2.0, 4.0, 6.0}; // Y = 2*X, dY/dX == 2 everywhere

  std::vector<real> dY_dx = compute_finite_difference(X, Y);

  for (real d : dY_dx)
  {
    CHECK(d == Approx(2.0));
  }
}
