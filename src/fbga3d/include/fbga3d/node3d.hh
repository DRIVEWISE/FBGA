#pragma once

#include <utils/types.hh>
#include <utils/gg_utils.hh>

namespace fb::fbga3d
{

using fb::utils::real;
using fb::utils::size_type;
using fb::utils::QUIET_NAN;
using fb::utils::SegmentType;

struct NodeStruct3D
{
  // STATIC MEMBERS (GIVEN)
  //// Length
  real s{0};
  //// geometry (given)
  ///// euler angles
  real mu{0.0};
  real phi{0.0};
  real theta{0.0};
  ///// euler derivatives
  real mu_prime{0.0};
  real phi_prime{0.0};
  real theta_prime{0.0};
  ///// euler second derivatives
  real mu_double_prime{0.0};
  real phi_double_prime{0.0};
  real theta_double_prime{0.0};
  ////// offset
  real n{0.0};
  real chi{0.0};
  ////// additional offsets
  real chi_prime{0.0};
  // STATIC MEMBERS (COMPUTED)
  //// Gravity corrections
  real g_x{0.0};
  real g_y{0.0};
  real g_z{0.0};
  //// Geometric Omegas
  real Omega_x{0.0};
  real Omega_y{0.0};
  real Omega_z{0.0};
  //// Geometric Omegas prime
  real Omega_x_prime{0.0};
  real Omega_y_prime{0.0};
  real Omega_z_prime{0.0};
  //
  real V_max{130.0};
  // Adherence scaling
  real alpha{1.0};
};

struct CellStruct3D
{
  real s_0{0};
  real s_1{0};
  real L{0};
  size_type ID0{0};
  size_type ID1{0};
  SegmentType m_type{SegmentType::UNKNOWN};
  real V_max{130.0};
  real V_dot{QUIET_NAN};
  real V_0{0.0};
  real V_1{0.0};
};

} // namespace fb::fbga3d
