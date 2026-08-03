#pragma once

#include <utils/types.hh>
#include <utils/gg_utils.hh>
#include <fbga3d/types3d.hh>

namespace fb::fbga3d
{

using fb::utils::real;
using fb::utils::BilinearInterpolator;
using fb::utils::LinearInterpolator;

class GggvIndy
{
private:
  BilinearInterpolator m_ay_max_bilinear;
  BilinearInterpolator m_ax_min_bilinear;
  BilinearInterpolator m_ax_max_bilinear;

  LinearInterpolator m_ax_eng_interpolator;

  ScalingGggvFactors m_scaling_factors;

public:
  GggvIndy();
  explicit GggvIndy(ScalingGggvFactors const &scaling_factors);
  explicit GggvIndy(SplineDataCollection const &spline_data);

  void set_scaling_factors(ScalingGggvFactors const &scaling_factors);
  void set_engine_max(EngineMaxGgv const &engine_max_ggv);

  // Combined (envelope + rho-blend + engine cap) longitudinal bounds. `alpha` scales the raw
  // GG envelope and the lateral limit used for the blend ratio -- NOT the engine cap, since
  // engine power doesn't depend on tire adherence. See FBGA3D_INTEGRATION_PLAN.md.
  [[nodiscard]] real a_x_push(real ay_tilde, real V, real az_tilde, real alpha = 1.0) const;
  [[nodiscard]] real a_x_pull(real ay_tilde, real V, real az_tilde, real alpha = 1.0) const;

  // Acceleration at the push/pull envelope split (zero net driver input). The spline data this
  // model is built from doesn't encode this split point, so it defaults to 0.
  [[nodiscard]] real a_x_neutral(real V) const;

  void setup_std();
  [[nodiscard]] real rho(real V, real az_tilde) const;
  [[nodiscard]] real rho_max(real V, real az_tilde) const;
  [[nodiscard]] real rho_min(real V, real az_tilde) const;

  [[nodiscard]] real a_x_eng(real V) const;
  [[nodiscard]] real a_x_max(real V, real az_tilde) const;
  [[nodiscard]] real a_x_min(real V, real az_tilde) const;
  [[nodiscard]] real a_y_lim(real V, real az_tilde) const;
};

} // namespace fb::fbga3d
