#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <cmath>

#include <solvers/brent_dekker.hh>

using namespace fb::solvers;
using Catch::Approx;

TEST_CASE("BrentDekker solves a simple linear function", "[solvers][BrentDekker]")
{
  BrentDekker solver(1e-12, 100);
  real x0 = 0.0;

  const bool ok = solver.solve([](real x) { return x - 3.0; }, 0.0, 10.0, x0);

  REQUIRE(ok);
  CHECK(x0 == Approx(3.0));
}

TEST_CASE("BrentDekker solves a quadratic function", "[solvers][BrentDekker]")
{
  BrentDekker solver(1e-12, 100);
  real x0 = 0.0;

  const bool ok = solver.solve([](real x) { return x * x - 2.0; }, 0.0, 2.0, x0);

  REQUIRE(ok);
  CHECK(x0 == Approx(std::sqrt(2.0)));
}

TEST_CASE("BrentDekker solves a transcendental function (Dottie number)", "[solvers][BrentDekker]")
{
  BrentDekker solver(1e-12, 200);
  real x0 = 0.0;

  const bool ok = solver.solve([](real x) { return std::cos(x) - x; }, 0.0, 1.0, x0);

  REQUIRE(ok);
  CHECK(x0 == Approx(0.7390851332151607));
}

TEST_CASE("BrentDekker returns the root immediately if it sits at a border", "[solvers][BrentDekker]")
{
  BrentDekker solver(1e-12, 100);
  real x0 = 0.0;

  const bool ok = solver.solve([](real x) { return x; }, 0.0, 1.0, x0);

  REQUIRE(ok);
  CHECK(x0 == Approx(0.0));
}

TEST_CASE("BrentDekker reports failure when the root is not bracketed", "[solvers][BrentDekker]")
{
  BrentDekker solver(1e-12, 100);
  real x0 = 0.0;

  // x^2 + 1 has no real root, so f(a) and f(b) share the same sign everywhere.
  const bool ok = solver.solve([](real x) { return x * x + 1.0; }, -1.0, 1.0, x0);

  CHECK_FALSE(ok);
}
