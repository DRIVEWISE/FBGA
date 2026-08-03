#pragma once

#include <utils/types.hh>

namespace fb::fbga3d
{

using fb::utils::real;

// Aerodynamic / geometric data for a motorcycle GG-diagram model
struct MotoData
{
  real b = 0.73;    // wheelbase
  real L_W = 1.5;   // length of the wheel
  real h = 0.69;    // height of the center of mass
  real mu_X = 1.30; // longitudinal friction coefficient
  real mu_Y = 1.40; // lateral friction coefficient
  real c_a_0 = 0.00; // aerodynamic coefficient
  real c_a_1 = 0.00; // aerodynamic coefficient
  real c_a_2 = 0.0 * 0.5 * 1.2 * 0.25 / (250.0 * 9.81); // aerodynamic coefficient
  real h_a = 0.51;  // height of the aerodynamic center
  real g = 9.81;    // gravitational acceleration
  real M = 250.0;   // mass of the motorcycle
  real P = 145.0 * 1000.0; // maximum power of the engine in Watts
};

class GggvMoto
{
private:
  MotoData m_aero_data; // Aerodynamic data

public:
  GggvMoto();

  // Combined (adherence/wheeling/engine-limited) longitudinal bounds. `alpha` scales the whole
  // result, including this model's own internal engine cap -- see FBGA3D_INTEGRATION_PLAN.md.
  [[nodiscard]] real a_x_push(real ay_tilde, real V, real az_tilde, real alpha = 1.0) const;
  [[nodiscard]] real a_x_pull(real ay_tilde, real V, real az_tilde, real alpha = 1.0) const;

  // Acceleration at the push/pull envelope split (zero net driver input): analytically this
  // coincides with the aero-drag deceleration for this closed-form model.
  [[nodiscard]] real a_x_neutral(real V) const;

  [[nodiscard]] real a_x_eng(real V) const;
  [[nodiscard]] real a_y_lim(real V, real az_tilde) const;

  [[nodiscard]] real a_x_aero(real V) const;

  void setup_std();

  [[nodiscard]] real constraint_ax_max_wheeling(real ay_hat, real az_hat, real v) const;
  [[nodiscard]] real constraint_ax_max_adherence(real ay_hat, real az_hat, real v) const;
  [[nodiscard]] real constraint_ax_min_stoppie(real ay_hat, real az_hat, real v) const;
  [[nodiscard]] real constraint_ax_min_adherence(real ay_hat, real az_hat, real v) const;
  [[nodiscard]] real constraint_ay_lim_adherence(real az_hat, real v) const;
};

} // namespace fb::fbga3d
