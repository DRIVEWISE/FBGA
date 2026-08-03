#include <fbga3d/gggv_moto.hh>

#include <algorithm>
#include <cmath>

namespace fb::fbga3d
{

void
GggvMoto::setup_std()
{
  // SETUP values
}

// -------------------------------------------------------------------

GggvMoto::GggvMoto()
{
  this->setup_std();
}

// -------------------------------------------------------------------

real
GggvMoto::a_x_push(real ay_tilde, real V, real az_tilde, real alpha) const
{
  const real a_y_lim = alpha * this->a_y_lim(V, az_tilde);
  const real ay_clip = std::clamp(ay_tilde, -a_y_lim, a_y_lim);
  real const a_x_max_adherence = this->constraint_ax_max_adherence(ay_clip, az_tilde, V);
  real const a_x_max_wheeling = this->constraint_ax_max_wheeling(ay_clip, az_tilde, V);
  real const a_x_max_engine = this->a_x_eng(V);
  return alpha * std::min({a_x_max_adherence, a_x_max_wheeling, a_x_max_engine});
}

real
GggvMoto::a_x_pull(real ay_tilde, real V, real az_tilde, real alpha) const
{
  const real a_y_lim = alpha * this->a_y_lim(V, az_tilde);
  const real ay_clip = std::clamp(ay_tilde, -a_y_lim, a_y_lim);
  real const a_x_min_adherence = this->constraint_ax_min_adherence(ay_clip, az_tilde, V);
  real const a_x_min_stoppie = this->constraint_ax_min_stoppie(ay_clip, az_tilde, V);
  return alpha * std::max(a_x_min_adherence, a_x_min_stoppie);
}

real
GggvMoto::a_x_neutral(real V) const
{
  return -this->a_x_aero(V);
}

real
GggvMoto::a_x_eng(real V) const
{
  const auto &M = this->m_aero_data.M;
  const auto &P = this->m_aero_data.P;
  const auto &c_a_0 = this->m_aero_data.c_a_0;
  const auto &c_a_1 = this->m_aero_data.c_a_1;
  const auto &c_a_2 = this->m_aero_data.c_a_2;
  const auto &g = this->m_aero_data.g;
  return P / (M * V) - (c_a_2 * V * V + c_a_1 * V + c_a_0) * g;
}

real
GggvMoto::a_x_aero(real V) const
{
  const auto &c_a_0 = this->m_aero_data.c_a_0;
  const auto &c_a_1 = this->m_aero_data.c_a_1;
  const auto &c_a_2 = this->m_aero_data.c_a_2;
  const auto &g = this->m_aero_data.g;
  const auto &M = this->m_aero_data.M;
  return -(c_a_2 * V * V + c_a_1 * V + c_a_0) * g * M;
}

real
GggvMoto::a_y_lim(real V, real az_tilde) const
{
  return this->constraint_ay_lim_adherence(az_tilde, V);
}

real
GggvMoto::constraint_ax_max_wheeling(real ay_hat, real az_hat, real v) const
{
  auto const &c_a_0 = this->m_aero_data.c_a_0;
  auto const &c_a_1 = this->m_aero_data.c_a_1;
  auto const &c_a_2 = this->m_aero_data.c_a_2;
  auto const &h = this->m_aero_data.h;
  auto const &h_a = this->m_aero_data.h_a;
  auto const &g = this->m_aero_data.g;
  auto const &b = this->m_aero_data.b;
  const real t1 = c_a_2 * g;
  const real t2 = v * v;
  const real t7 = c_a_1 * g;
  const real t12 = c_a_0 * g;
  const real t15 = ay_hat * ay_hat;
  const real t16 = az_hat * az_hat;
  const real t18 = sqrt(t15 + t16);
  return (0.1e1 / h * (-h * t1 * t2 - h * t7 * v + h_a * t1 * t2 + h_a * t7 * v + b * t18 - h * t12 + h_a * t12));
}

real
GggvMoto::constraint_ax_max_adherence(real ay_hat, real az_hat, real v) const
{
  auto const &c_a_0 = this->m_aero_data.c_a_0;
  auto const &c_a_1 = this->m_aero_data.c_a_1;
  auto const &c_a_2 = this->m_aero_data.c_a_2;
  auto const &h = this->m_aero_data.h;
  auto const &h_a = this->m_aero_data.h_a;
  auto const &g = this->m_aero_data.g;
  auto const &b = this->m_aero_data.b;
  auto const &L_W = this->m_aero_data.L_W;
  auto const &mu_X = this->m_aero_data.mu_X;
  auto const &mu_Y = this->m_aero_data.mu_Y;
  const real t1 = mu_X * mu_X;
  const real t2 = t1 * h;
  const real t3 = az_hat * mu_Y;
  const real t4 = t3 - ay_hat;
  const real t6 = t3 + ay_hat;
  const real t7 = b - L_W;
  const real t9 = ay_hat * ay_hat;
  const real t10 = az_hat * az_hat;
  const real t11 = t9 + t10;
  const real t12 = sqrt(t11);
  const real t15 = L_W * L_W;
  const real t16 = mu_Y * mu_Y;
  const real t17 = t15 * t16;
  const real t19 = (h - h_a) * t2;
  const real t22 = v * v;
  const real t25 = c_a_1 * v + c_a_2 * t22 + c_a_0;
  const real t41 = t7 * t7;
  const real t44 = g * g;
  const real t45 = h_a * h_a;
  const real t47 = t25 * t25;
  const real t53 = sqrt(
    t1 * (0.2e1 * g * h_a * t12 * t25 * t7 + t44 * t45 * t47 + t10 * t41 + t41 * t9) * t15 * t4 * t16 * t11 * t6
  );
  const real t55 = h * h;
  const real t56 = t1 * t55;
  return (
    0.1e1 / (t9 * (-t17 - t56) + (t56 - t15) * t16 * t10) *
    (t12 * t7 * t6 * t4 * t2 + t9 * t25 * (t17 + t19) * g - t10 * t16 * t25 * g * (-t15 + t19) - t53)
  );
}

real
GggvMoto::constraint_ax_min_stoppie(real ay_hat, real az_hat, real v) const
{
  auto const &c_a_0 = this->m_aero_data.c_a_0;
  auto const &c_a_1 = this->m_aero_data.c_a_1;
  auto const &c_a_2 = this->m_aero_data.c_a_2;
  auto const &h = this->m_aero_data.h;
  auto const &h_a = this->m_aero_data.h_a;
  auto const &g = this->m_aero_data.g;
  auto const &b = this->m_aero_data.b;
  auto const &L_W = this->m_aero_data.L_W;
  real const t1 = c_a_2 * g;
  real const t2 = v * v;
  real const t7 = c_a_1 * g;
  real const t12 = c_a_0 * g;
  real const t15 = ay_hat * ay_hat;
  real const t16 = az_hat * az_hat;
  real const t18 = sqrt(t15 + t16);
  return (-0.1e1 / h * (h * t1 * t2 + h * t7 * v - h_a * t1 * t2 - h_a * t7 * v + L_W * t18 - b * t18 + h * t12 - h_a * t12));
}

real
GggvMoto::constraint_ax_min_adherence(real ay_hat, real az_hat, real v) const
{
  auto const &c_a_0 = this->m_aero_data.c_a_0;
  auto const &c_a_1 = this->m_aero_data.c_a_1;
  auto const &c_a_2 = this->m_aero_data.c_a_2;
  auto const &g = this->m_aero_data.g;
  auto const &mu_X = this->m_aero_data.mu_X;
  auto const &mu_Y = this->m_aero_data.mu_Y;
  real const t2 = v * v;
  real const t10 = mu_Y * mu_Y;
  real const t11 = mu_X * mu_X;
  real const t13 = az_hat * az_hat;
  real const t15 = ay_hat * ay_hat;
  real const t18 = sqrt(t10 * t11 * t13 - t15 * t11);
  return (0.1e1 / mu_Y * (-c_a_1 * g * mu_Y * v - c_a_2 * g * mu_Y * t2 - c_a_0 * g * mu_Y - t18));
}

// --------------------------------------------------------------------------------------------

real
GggvMoto::constraint_ay_lim_adherence(real az_hat, real /*v*/) const
{
  auto const &mu_Y = this->m_aero_data.mu_Y;
  return az_hat * mu_Y;
}

} // namespace fb::fbga3d
