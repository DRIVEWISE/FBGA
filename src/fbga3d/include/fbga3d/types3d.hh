#pragma once

#include <utils/types.hh>

#include <string>
#include <vector>

namespace fb::fbga3d
{

using fb::utils::real;
using fb::utils::floating;
using fb::utils::integer;
using fb::utils::GRAVITY;
using fb::utils::INFTY;
using fb::utils::QUIET_NAN;

constexpr real VMAX = 130.0; // Maximum velocity

struct SolverParams
{
  real tolerance        = fb::utils::EPSILON_HIGH * 10; // matches FBGA_3D's STD_TOL (1e-10)
  int max_iter           = 200;                          // matches FBGA_3D's STD_MAX_ITER
  std::string verbosity  = "zero";                       // matches FBGA_3D's STD_VERBOSE
};

struct ScalingGggvFactors
{
  real ax_max_scale = 1.0;       // Scale factor for a_x_max
  real ax_min_scale = 1.0;       // Scale factor for a_x_min
  real ay_scale = 1.0;           // Scale factor for a_y
  real gg_exponent_ax_pos = 1.3; // Exponent for a_x_max positive
  real gg_exponent_ax_neg = 1.3; // Exponent for a_x_max negative
};

struct SplineDataCollection
{
  std::vector<real> v_data{0.0, 90.0};
  std::vector<real> az_data{5.0, 15.0};
  std::vector<real> ay_max_data{10.0, 10.0, 10.0, 10.0};
  std::vector<real> ax_min_data{-10.0, -10.0, -10.0, -10.0};
  std::vector<real> ax_max_data{10.0, 10.0, 10.0, 10.0};
  ScalingGggvFactors scaling_factors;
  std::vector<real> v_data_eng{0.0, 90.0};   // Velocity data for engine max GGV
  std::vector<real> ax_eng_data{13.0, 0.0};  // Acceleration data for engine max GGV
};

struct EngineMaxGgv
{
  std::vector<real> v_data_eng{0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0};
  std::vector<real> ax_eng_data{5.0, 10.0, 13.0, 10.0, 7.0, 4.6, 2.5, 0.3, 0.0, 0.0};
};

// Road angles and their geometric derivatives, sampled along the path
struct RoadAnglesAndDerivativesContainer
{
  std::vector<real> mu;           // Pitch angle
  std::vector<real> phi;          // Roll angle
  std::vector<real> theta;        // Yaw angle
  std::vector<real> mu_prime;     // Geometric derivative of the pitch angle
  std::vector<real> phi_prime;    // Geometric derivative of the roll angle
  std::vector<real> theta_prime;  // Geometric derivative of the yaw angle
  std::vector<real> abscissa;     // Arc-length abscissa
};

struct TrajectoryOffsetContainer
{
  std::vector<real> n;
  std::vector<real> chi;
};

struct AdherenceContainer
{
  std::vector<real> alpha; // Adherence coefficient
};

struct TrajectoryOffsetAndAnglesContainer
{
  TrajectoryOffsetContainer offset;
  RoadAnglesAndDerivativesContainer reference;
  AdherenceContainer adherence;
};

struct SolutionContainer
{
  std::vector<real> V0;    // Initial velocity
  std::vector<real> V1;    // Final velocity
  std::vector<real> V_dot; // Velocity derivative
};

struct OutputPlotGgvShell
{
  std::vector<floating> a_tilde_x_max;
  std::vector<floating> a_tilde_x_min;
  std::vector<floating> a_tilde_y;
  std::vector<floating> v;
};

struct InputPlotGgvShell
{
  std::vector<real> v = {0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0};
  real az = GRAVITY;
  real alpha = 1.0;
};

struct OutputPlotGgvTraj
{
  std::vector<floating> a_tilde_x;
  std::vector<floating> a_tilde_y;
  std::vector<floating> v;
};

struct YellowFlagData
{
  real v_des_max = VMAX;
  real a_des_min = -100.0;
  bool is_yellow = false; // Flag to indicate if the yellow flag is active
};

struct AblationFlags
{
  bool compute_chi_dot = true;
  bool compute_s_ddot = true;
  bool compute_Omega_xyz_prime = true;
  bool compute_Euler_prime = true;
  bool compute_Euler_double_prime = true;
};

struct ConstraintViolation
{
  real maximum = -INFTY;
  real minimum = +INFTY;
  real average = QUIET_NAN;
  integer id_max = -1;
  integer id_min = -1;
  real H_inf = INFTY;
  real cumulative = 0.0;
  std::vector<real> values;
  std::vector<real> values_violation;
  integer num_violations = 0;
};

struct ConstraintViolationForSegments
{
  std::vector<std::vector<real>> values;
  std::vector<std::vector<real>> values_violation;
  real maximum = -INFTY;
  real minimum = INFTY;
  real average = QUIET_NAN;
  integer id_max = -1;
  integer id_min = -1;
  real H_inf = INFTY;
  real cumulative = 0.0;
  std::vector<real> cumulative_by_segment;
  std::vector<real> average_by_segment;
  std::vector<real> H_inf_by_segment;
  integer num_violations = 0;
};

// Per-node kinematic quantities computed during a single evaluation step
struct Context
{
  real cos_chi{1}, sin_chi{0};
  real den_common{1}, inv_den_common{1}, inv_den_common2{1}, inv_den_common3{1};
  real s_dot{0}, w{0}, chi_dot{0}, n_dot{0};

  real s_ddotA{0}, s_ddotB{0}, s_ddot{0}, w_dot{0};

  real omega_hat_x{0}, omega_hat_y{0}, omega_hat_z{0};

  real a_tilde_x{0}, a_tilde_y{0}, a_tilde_z{0};

  real a_tilde_y_lim{0}, a_tilde_y_clip{0};
  // Already-combined (envelope + engine cap, per-model alpha convention baked in via
  // GggvModel::a_x_push/a_x_pull) longitudinal bounds -- see FBGA3D_INTEGRATION_PLAN.md.
  real a_tilde_x_max_gg{0}, a_tilde_x_min_gg{0}, a_tilde_x_eng{0};

  real a_tilde_x_max{0}, a_tilde_x_min{0};

  bool at_lateral_limit{false};
};

} // namespace fb::fbga3d
