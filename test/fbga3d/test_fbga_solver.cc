#include <catch2/catch_test_macros.hpp>

#include <fbga3d/fbga_indy.hh>
#include <fbga3d/fbga_moto.hh>

#include <cmath>

using namespace fb::fbga3d;

namespace
{

// Flat, non-banked track (mu = phi = 0 everywhere) with a simple curvature profile via
// theta_prime, mirroring examples/fbga/basic.cc's SS/KK arrays for the 2D solver.
TrajectoryOffsetAndAnglesContainer make_flat_trajectory()
{
  TrajectoryOffsetAndAnglesContainer TOA;
  TOA.reference.abscissa = {0.0, 50.0, 100.0, 150.0, 200.0};
  TOA.reference.mu = {0.0, 0.0, 0.0, 0.0, 0.0};
  TOA.reference.phi = {0.0, 0.0, 0.0, 0.0, 0.0};
  TOA.reference.theta = {0.0, 0.05, 0.55, 0.8, 0.8};
  TOA.reference.mu_prime = {0.0, 0.0, 0.0, 0.0, 0.0};
  TOA.reference.phi_prime = {0.0, 0.0, 0.0, 0.0, 0.0};
  TOA.reference.theta_prime = {0.001, 0.01, 0.005, -0.01, 0.0};
  TOA.offset.n = {0.0, 0.0, 0.0, 0.0, 0.0};
  TOA.offset.chi = {0.0, 0.0, 0.0, 0.0, 0.0};
  TOA.adherence.alpha = {1.0, 1.0, 1.0, 1.0, 1.0};
  return TOA;
}

} // namespace

TEST_CASE("FbgaMoto::compute returns a finite total time on a flat track", "[fbga3d][fbga_moto]")
{
  FbgaMoto solver;
  const auto TOA = make_flat_trajectory();
  const real total_time = solver.compute(TOA, 20.0);

  CHECK(std::isfinite(total_time));
  CHECK(total_time > 0.0);

  // A cell can legitimately fail to converge (e.g. a curvature/speed combination outside
  // GggvMoto's adherence envelope) -- the solver falls back to a safe clamped value rather
  // than propagating NaN. What matters is that the fallback keeps the result finite.
  const auto sol = solver.get_solution();
  CHECK(sol.V0.size() == TOA.reference.abscissa.size() - 1);
  for (real v : sol.V0) CHECK(std::isfinite(v));
  for (real v : sol.V1) CHECK(std::isfinite(v));
}

TEST_CASE("FbgaMoto eval_V/eval_A_tilde_x accessors are finite along the solved profile", "[fbga3d][fbga_moto]")
{
  FbgaMoto solver;
  const auto TOA = make_flat_trajectory();
  solver.compute(TOA, 20.0);

  for (real s = 1.0; s < 199.0; s += 20.0)
  {
    CHECK(std::isfinite(solver.eval_V(s)));
    CHECK(std::isfinite(solver.eval_A_tilde_x(s)));
    CHECK(std::isfinite(solver.eval_Vmax(s)));
  }
}

// GggvIndy::setup_std() (invoked by FbgaIndy's default-constructed model) reads
// "./data/INDY/*.npy" relative to the working directory -- see this target's
// WORKING_DIRECTORY in test/fbga3d/CMakeLists.txt.
TEST_CASE("FbgaIndy::compute returns a finite total time on a flat track", "[fbga3d][fbga_indy]")
{
  FbgaIndy solver;
  const auto TOA = make_flat_trajectory();
  const real total_time = solver.compute(TOA, 20.0);

  CHECK(std::isfinite(total_time));
  CHECK(total_time > 0.0);

  const auto sol = solver.get_solution();
  CHECK(sol.V0.size() == TOA.reference.abscissa.size() - 1);
  for (real v : sol.V0) CHECK(std::isfinite(v));
  for (real v : sol.V1) CHECK(std::isfinite(v));
}
