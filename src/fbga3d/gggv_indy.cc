#include <fbga3d/gggv_indy.hh>

#include <npy.hpp>

#include <algorithm>
#include <cmath>
#include <string>
#include <utility>
#include <vector>

#define DEBUG_SPLINE 0

namespace fb::fbga3d
{

void
GggvIndy::setup_std()
{
  const std::string base_folder = "./data/INDY/";
  const std::string g_list_path = base_folder + "g_list.npy";
  const std::string v_list_path = base_folder + "v_list.npy";
  const std::string ay_max_path = base_folder + "ay_max.npy";
  const std::string ax_min_path = base_folder + "ax_min.npy";
  const std::string ax_max_path = base_folder + "ax_max.npy";

  const npy::npy_data g_npy_data = npy::read_npy<real>(g_list_path);
  const npy::npy_data v_npy_data = npy::read_npy<real>(v_list_path);
  const npy::npy_data ay_max_npy_data = npy::read_npy<real>(ay_max_path);
  const npy::npy_data ax_min_npy_data = npy::read_npy<real>(ax_min_path);
  const npy::npy_data ax_max_npy_data = npy::read_npy<real>(ax_max_path);

  const std::vector<real> g_data = std::move(g_npy_data.data);
  const std::vector<real> v_data = std::move(v_npy_data.data);
  const std::vector<real> ay_max_data = std::move(ay_max_npy_data.data);
  const std::vector<real> ax_min_data = std::move(ax_min_npy_data.data);
  const std::vector<real> ax_max_data = std::move(ax_max_npy_data.data);

#if DEBUG_SPLINE
  std::vector<unsigned long> g_shape = g_npy_data.shape;
  std::vector<unsigned long> v_shape = v_npy_data.shape;
  std::vector<unsigned long> ay_max_shape = ay_max_npy_data.shape;
  std::vector<unsigned long> ax_min_shape = ax_min_npy_data.shape;
  std::vector<unsigned long> ax_max_shape = ax_max_npy_data.shape;
  bool g_fortran_order = g_npy_data.fortran_order;
  bool v_fortran_order = v_npy_data.fortran_order;
  bool ay_max_fortran_order = ay_max_npy_data.fortran_order;
  bool ax_min_fortran_order = ax_min_npy_data.fortran_order;
  bool ax_max_fortran_order = ax_max_npy_data.fortran_order;
#endif

  this->m_ay_max_bilinear.setup(v_data, g_data, ay_max_data);
  this->m_ax_min_bilinear.setup(v_data, g_data, ax_min_data);
  this->m_ax_max_bilinear.setup(v_data, g_data, ax_max_data);

  EngineMaxGgv eng_ggv;
  this->m_ax_eng_interpolator.setup(eng_ggv.v_data_eng, eng_ggv.ax_eng_data);
}

// -------------------------------------------------------------------

void
GggvIndy::set_scaling_factors(ScalingGggvFactors const &scaling_factors)
{
  this->m_scaling_factors = scaling_factors;
}

// -------------------------------------------------------------------

void
GggvIndy::set_engine_max(EngineMaxGgv const &engine_max_ggv)
{
  this->m_ax_eng_interpolator.setup(engine_max_ggv.v_data_eng, engine_max_ggv.ax_eng_data);
}

// -------------------------------------------------------------------

GggvIndy::GggvIndy()
{
  this->setup_std();
}

// -------------------------------------------------------------------

GggvIndy::GggvIndy(ScalingGggvFactors const &scaling_factors)
{
  this->setup_std();
  this->set_scaling_factors(scaling_factors);
}

// -------------------------------------------------------------------

GggvIndy::GggvIndy(SplineDataCollection const &spline_data)
{
  this->m_ay_max_bilinear.setup(spline_data.v_data, spline_data.az_data, spline_data.ay_max_data);
  this->m_ax_min_bilinear.setup(spline_data.v_data, spline_data.az_data, spline_data.ax_min_data);
  this->m_ax_max_bilinear.setup(spline_data.v_data, spline_data.az_data, spline_data.ax_max_data);

  this->m_ax_eng_interpolator.setup(spline_data.v_data_eng, spline_data.ax_eng_data);

  this->set_scaling_factors(spline_data.scaling_factors);
}

// -------------------------------------------------------------------

real
GggvIndy::a_x_push(real ay_tilde, real V, real az_tilde, real alpha) const
{
  const real a_x_max_raw = this->a_x_max(V, az_tilde);
  const real a_y_lim = alpha * this->a_y_lim(V, az_tilde);
  const real rho_max = this->rho_max(V, az_tilde);
  const real a_x_eng = this->a_x_eng(V);
  const real ay_tilde_abs = std::abs(ay_tilde);
  const real ay_tilde_absclip = std::min(ay_tilde_abs, a_y_lim);
  // Matches SolverParams::tolerance's default (1e-10), not the 1e-12 gggv_INDY::a_x_push used in
  // FBGA_3D -- that version was never actually called by FBGA_INDY, which used this wider margin
  // in its own inline blend instead. See FBGA3D_INTEGRATION_PLAN.md.
  if (std::abs(ay_tilde_abs - a_y_lim) < 1e-10)
  {
    return std::min(0.0, a_x_eng);
  }
  return std::min(
    std::pow((1 - std::pow(ay_tilde_absclip / a_y_lim, rho_max)), (1 / rho_max)) * alpha * a_x_max_raw,
    a_x_eng
  );
}

real
GggvIndy::a_x_pull(real ay_tilde, real V, real az_tilde, real alpha) const
{
  const real a_x_min_raw = this->a_x_min(V, az_tilde);
  const real a_y_lim = alpha * this->a_y_lim(V, az_tilde);
  const real rho_min = this->rho_min(V, az_tilde);
  const real ay_tilde_abs = std::abs(ay_tilde);
  const real ay_tilde_absclip = std::min(ay_tilde_abs, a_y_lim);

  // Matches SolverParams::tolerance's default (1e-10), not the 1e-12 gggv_INDY::a_x_push used in
  // FBGA_3D -- that version was never actually called by FBGA_INDY, which used this wider margin
  // in its own inline blend instead. See FBGA3D_INTEGRATION_PLAN.md.
  if (std::abs(ay_tilde_abs - a_y_lim) < 1e-10)
  {
    return std::max(0.0, alpha * a_x_min_raw);
  }
  return std::pow((1 - std::pow(ay_tilde_absclip / a_y_lim, rho_min)), (1 / rho_min)) * alpha * a_x_min_raw;
}

real
GggvIndy::a_x_neutral(real /*V*/) const
{
  return 0.0;
}

real
GggvIndy::rho(real /*V*/, real /*az_tilde*/) const
{
  return 1.3;
}

real
GggvIndy::rho_max(real /*V*/, real /*az_tilde*/) const
{
  return this->m_scaling_factors.gg_exponent_ax_pos;
}

real
GggvIndy::rho_min(real /*V*/, real /*az_tilde*/) const
{
  return this->m_scaling_factors.gg_exponent_ax_neg;
}

real
GggvIndy::a_x_eng(real V) const
{
  return this->m_ax_eng_interpolator.eval(V);
}

real
GggvIndy::a_x_max(real V, real az_tilde) const
{
  return this->m_ax_max_bilinear.eval(V, az_tilde) * this->m_scaling_factors.ax_max_scale;
}

real
GggvIndy::a_x_min(real V, real az_tilde) const
{
  return this->m_ax_min_bilinear.eval(V, az_tilde) * this->m_scaling_factors.ax_min_scale;
}

real
GggvIndy::a_y_lim(real V, real az_tilde) const
{
  return this->m_ay_max_bilinear.eval(V, az_tilde) * this->m_scaling_factors.ay_scale;
}

} // namespace fb::fbga3d
