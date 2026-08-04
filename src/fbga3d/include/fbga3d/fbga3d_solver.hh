#pragma once

#include <utils/gg_utils.hh>
#include <utils/types.hh>
#include <solvers/brent_dekker.hh>
#include <fbga3d/node3d.hh>
#include <fbga3d/types3d.hh>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

namespace fb::fbga3d
{

using fb::utils::clip;
using fb::utils::GRAVITY;
using fb::utils::integer;
using fb::utils::QUIET_NAN;
using fb::utils::real;
using fb::utils::SegmentType;
using fb::utils::size_type;
using fb::solvers::BrentDekker;

// Per-vehicle default tolerances. Specialize for a given GggvModel *before* instantiating
// Fbga3dSolver<Model> to override -- see fbga_moto.hh for FBGA_3D's FBGA_MOTO defaults.
template <class GggvModel>
struct Fbga3dDefaults
{
  static constexpr real lat_tol = 0.0;
  static constexpr real lat_tol_vmax = 0.0;
};

// Forward-backward 3D velocity-profile solver, generic over the vehicle's GG-diagram model
// (GggvIndy, GggvMoto, ...). Ported from FBGA_3D's FBGA_INDY/FBGA_MOTO, which were near-identical
// hand-duplicated classes; see FBGA3D_INTEGRATION_PLAN.md for the reconciliation of the handful
// of places they actually differed (alpha handling, aero-neutral point, lateral tolerances).
//
// Core solve path only: plotting helpers (eval_shell_plot/eval_shell_plotpy) and the
// constraint-satisfaction-checking machinery from the original classes are not ported.
template <class GggvModel>
class Fbga3dSolver
{
public:
  static constexpr size_type DEFAULT_RESERVE = 100;

  Fbga3dSolver()
  {
    this->Nodes.reserve(DEFAULT_RESERVE);
    this->Cells.reserve(DEFAULT_RESERVE);
    this->dump_seg_id.reserve(DEFAULT_RESERVE);
  }

  [[nodiscard]] GggvModel &model() { return this->gggv; }
  [[nodiscard]] GggvModel const &model() const { return this->gggv; }

  void set_max_velocity(real velocity) { this->max_velocity = velocity; }

  // --------------------------------------------------------------------------------------------

  real compute(
    TrajectoryOffsetAndAnglesContainer const &TOA,
    real V0 = -1.0,
    YellowFlagData const &yellow_flag = YellowFlagData()
  )
  {
    this->m_yellow_flag_data = yellow_flag;
    this->dump_seg_id.clear();
    this->v_I = V0;
    if (this->v_I <= 0.0)
    {
      this->v_I = this->max_velocity;
    }
    this->create_nodes_cells(TOA);
    this->compute_Vmax();
    this->yellow_index = 0;
    if (this->m_yellow_flag_data.is_yellow)
    {
      this->BY();
    }
    this->FW();
    this->BW();
    return this->compute_total_time();
  }

  [[nodiscard]] SolutionContainer get_solution() const
  {
    SolutionContainer sol;
    const size_type N = this->Cells.size();
    sol.V0.resize(N);
    sol.V1.resize(N);
    sol.V_dot.resize(N);
    for (size_type i = 0; i < N; i++)
    {
      sol.V0[i] = this->Cells[i].V_0;
      sol.V1[i] = this->Cells[i].V_1;
      sol.V_dot[i] = this->Cells[i].V_dot;
    }
    return sol;
  }

  [[nodiscard]] std::vector<int> get_dump() const { return this->dump_seg_id_int(); }

  // --------------------------------------------------------------------------------------------
  // Accessors (by arc-length s / time t)

  [[nodiscard]] integer get_cell_idx(const real s) const
  {
    integer left = 0;
    integer right = static_cast<integer>(this->Cells.size()) - 1;
    if (s >= this->Cells[right].s_1)
    {
      return right;
    }
    if (s < this->Cells[left].s_0)
    {
      return -1;
    }
    while (left <= right)
    {
      const integer mid = left + (right - left) / 2;
      const auto &cell = this->Cells[mid];
      if (s >= cell.s_0 && s < cell.s_1)
      {
        return mid;
      }
      if (s < cell.s_0)
      {
        right = mid - 1;
      }
      else
      {
        left = mid + 1;
      }
    }
    return -1;
  }

  [[nodiscard]] real eval_Vmax(real s) const
  {
    const integer cell_idx = this->get_cell_idx(s);
    const auto &cell = this->Cells[cell_idx];
    const real s_norm = (s - cell.s_0) / (cell.s_1 - cell.s_0);
    const real V_max_0 = this->Nodes[cell.ID0].V_max;
    const real V_max_1 = this->Nodes[cell.ID1].V_max;
    return (V_max_0 * (1 - s_norm)) + (V_max_1 * s_norm);
  }

  [[nodiscard]] real eval_Omega_x(real s) const { return this->interp_node_field(s, &NodeStruct3D::Omega_x); }
  [[nodiscard]] real eval_Omega_y(real s) const { return this->interp_node_field(s, &NodeStruct3D::Omega_y); }
  [[nodiscard]] real eval_Omega_z(real s) const { return this->interp_node_field(s, &NodeStruct3D::Omega_z); }

  [[nodiscard]] real eval_V(real s) const
  {
    const integer cell_idx = this->get_cell_idx(s);
    const auto &cell = this->Cells[cell_idx];
    const auto &node = this->Nodes[cell.ID0];
    return eval_V(node, cell, s - cell.s_0);
  }

  [[nodiscard]] real eval_V_dot(real s) const
  {
    const integer cell_idx = this->get_cell_idx(s);
    return this->Cells[cell_idx].V_dot;
  }

  [[nodiscard]] real eval_A_hat_x(real s) const
  {
    const integer cell_idx = this->get_cell_idx(s);
    const auto &cell = this->Cells[cell_idx];
    const auto &node = this->Nodes[cell.ID0];
    const real V = eval_V(node, cell, s - cell.s_0);
    return this->eval_a_hat_x(node, V, cell.V_dot);
  }

  [[nodiscard]] real eval_A_hat_y(real s) const
  {
    const integer cell_idx = this->get_cell_idx(s);
    const auto &cell = this->Cells[cell_idx];
    const auto &node = this->Nodes[cell.ID0];
    const real V = eval_V(node, cell, s - cell.s_0);
    return this->eval_a_hat_y(node, V);
  }

  [[nodiscard]] real eval_A_hat_z(real s) const
  {
    const integer cell_idx = this->get_cell_idx(s);
    const auto &cell = this->Cells[cell_idx];
    const auto &node = this->Nodes[cell.ID0];
    const real V = eval_V(node, cell, s - cell.s_0);
    return this->eval_a_hat_z(node, V, cell.V_dot);
  }

  [[nodiscard]] real eval_A_tilde_x(real s) const
  {
    const integer cell_idx = this->get_cell_idx(s);
    const auto &cell = this->Cells[cell_idx];
    const auto &node = this->Nodes[cell.ID0];
    const real V = eval_V(node, cell, s - cell.s_0);
    return this->eval_a_tilde_x(node, V, cell.V_dot);
  }

  [[nodiscard]] real eval_A_tilde_y(real s) const
  {
    const integer cell_idx = this->get_cell_idx(s);
    const auto &cell = this->Cells[cell_idx];
    const auto &node = this->Nodes[cell.ID0];
    const real V = eval_V(node, cell, s - cell.s_0);
    return this->eval_a_tilde_y(node, V);
  }

  [[nodiscard]] real eval_A_tilde_z(real s) const
  {
    const integer cell_idx = this->get_cell_idx(s);
    const auto &cell = this->Cells[cell_idx];
    const auto &node = this->Nodes[cell.ID0];
    const real V = eval_V(node, cell, s - cell.s_0);
    return this->eval_a_tilde_z(node, V, cell.V_dot);
  }

  [[nodiscard]] real eval_g_x(real s) const { return this->interp_cell_field(s, &NodeStruct3D::g_x); }
  [[nodiscard]] real eval_g_y(real s) const { return this->interp_cell_field(s, &NodeStruct3D::g_y); }
  [[nodiscard]] real eval_g_z(real s) const { return this->interp_cell_field(s, &NodeStruct3D::g_z); }

  [[nodiscard]] SegmentType eval_segment_type(real s) const
  {
    return this->Cells[this->get_cell_idx(s)].m_type;
  }

  [[nodiscard]] real eval_alpha(real s) const
  {
    const integer cell_idx = this->get_cell_idx(s);
    const auto &cell = this->Cells[cell_idx];
    const auto &node_0 = this->Nodes[cell.ID0];
    const auto &node_1 = this->Nodes[cell.ID1];
    const real s_norm = (s - cell.s_0) / cell.L;
    return (node_0.alpha * (1 - s_norm)) + (node_1.alpha * s_norm);
  }

private:
  GggvModel gggv;
  BrentDekker BD;
  SolverParams solver_p;
  std::vector<NodeStruct3D> Nodes;
  std::vector<CellStruct3D> Cells;
  real v_I{0.0};
  real max_velocity{VMAX};
  std::vector<size_type> dump_seg_id;
  real m_lat_tol{Fbga3dDefaults<GggvModel>::lat_tol};
  real m_lat_tol_vmax{Fbga3dDefaults<GggvModel>::lat_tol_vmax};
  AblationFlags m_ablation_flags;
  real m_lateral_shrink_factor{1.0};
  integer yellow_index{0};
  YellowFlagData m_yellow_flag_data;

  [[nodiscard]] std::vector<int> dump_seg_id_int() const
  {
    return std::vector<int>(this->dump_seg_id.begin(), this->dump_seg_id.end());
  }

  template <class Member>
  [[nodiscard]] real interp_node_field(real s, Member NodeStruct3D::*field) const
  {
    const integer cell_idx = this->get_cell_idx(s);
    const auto &cell = this->Cells[cell_idx];
    const real s_norm = (s - cell.s_0) / (cell.s_1 - cell.s_0);
    const real v0 = this->Nodes[cell.ID0].*field;
    const real v1 = this->Nodes[cell.ID1].*field;
    return (v0 * (1 - s_norm)) + (v1 * s_norm);
  }

  template <class Member>
  [[nodiscard]] real interp_cell_field(real s, Member NodeStruct3D::*field) const
  {
    const integer cell_idx = this->get_cell_idx(s);
    const auto &cell = this->Cells[cell_idx];
    const real s_norm = (s - cell.s_0) / cell.L;
    const real v0 = this->Nodes[cell.ID0].*field;
    const real v1 = this->Nodes[cell.ID1].*field;
    return (v0 * (1 - s_norm)) + (v1 * s_norm);
  }

  // Everything derivable from (node, V) alone -- independent of V_dot. Root-finds that hold
  // (node, V) fixed while searching over V_dot (get_max/min_longitudinal_performance(_YF)) compute
  // this once and reuse it across every Brent-Dekker iteration, instead of recomputing the same
  // trig/kinematics on every evaluation. See FBGA3D_INTEGRATION_PLAN.md.
  struct NodeVelocityState
  {
    real cos_chi{1}, sin_chi{0};
    real inv_den_common{1}, inv_den_common2{1};
    real s_dot{0}, w{0}, chi_dot{0}, n_dot{0};
    real omega_hat_x{0}, omega_hat_y{0}, omega_hat_z{0};
    real a_tilde_y{0}; // V_dot-independent, so computed here once
    real s_ddotA{0};   // the V_dot-independent half of s_ddot
  };

  [[nodiscard]] NodeVelocityState compute_velocity_state(NodeStruct3D const &node, real V) const
  {
    NodeVelocityState st;
    st.cos_chi = std::cos(node.chi);
    st.sin_chi = std::sin(node.chi);
    const real den_common = node.Omega_z * node.n - 1.0;
    st.inv_den_common = 1.0 / den_common;
    st.inv_den_common2 = st.inv_den_common * st.inv_den_common;
    st.s_dot = V * st.cos_chi * (-st.inv_den_common);
    st.w = node.Omega_x * st.s_dot * node.n;
    st.chi_dot = this->m_ablation_flags.compute_chi_dot ? node.chi_prime * st.s_dot : 0.0;
    st.n_dot = V * st.sin_chi;
    st.omega_hat_x = (st.sin_chi * node.Omega_y + st.cos_chi * node.Omega_x) * st.s_dot;
    st.omega_hat_y = (st.cos_chi * node.Omega_y - st.sin_chi * node.Omega_x) * st.s_dot;
    st.omega_hat_z = node.Omega_z * st.s_dot + st.chi_dot;
    st.a_tilde_y = (st.omega_hat_z * V) - (st.omega_hat_x * st.w) + node.g_y;
    st.s_ddotA = this->m_ablation_flags.compute_s_ddot
      ? -V * st.cos_chi * (node.Omega_z_prime * st.s_dot * node.n + node.Omega_z * st.sin_chi * V) * st.inv_den_common2
      : 0.0;
    return st;
  }

  // --------------------------------------------------------------------------------------------

  void create_nodes_cells(TrajectoryOffsetAndAnglesContainer const &TOA)
  {
    this->Nodes.clear();
    this->Cells.clear();
    const size_type N = TOA.reference.abscissa.size();
    this->Nodes.resize(N);
    this->Cells.resize(N - 1);

    std::vector<real> tmp_dchi = fb::utils::compute_finite_difference(TOA.reference.abscissa, TOA.offset.chi);
    std::vector<real> tmp_dmu_prime =
      fb::utils::compute_finite_difference(TOA.reference.abscissa, TOA.reference.mu_prime);
    std::vector<real> tmp_dphi_prime =
      fb::utils::compute_finite_difference(TOA.reference.abscissa, TOA.reference.phi_prime);
    std::vector<real> tmp_dtheta_prime =
      fb::utils::compute_finite_difference(TOA.reference.abscissa, TOA.reference.theta_prime);

    for (size_type i = 0; i < N; i++)
    {
      auto &node = this->Nodes[i];
      node.s = TOA.reference.abscissa[i];
      node.mu = TOA.reference.mu[i];
      node.phi = TOA.reference.phi[i];
      node.theta = TOA.reference.theta[i];
      if (this->m_ablation_flags.compute_Euler_prime)
      {
        node.mu_prime = TOA.reference.mu_prime[i];
        node.phi_prime = TOA.reference.phi_prime[i];
        node.theta_prime = TOA.reference.theta_prime[i];
      }
      if (this->m_ablation_flags.compute_Euler_double_prime)
      {
        node.mu_double_prime = tmp_dmu_prime[i];
        node.phi_double_prime = tmp_dphi_prime[i];
        node.theta_double_prime = tmp_dtheta_prime[i];
      }
      node.n = TOA.offset.n[i];
      node.chi = TOA.offset.chi[i];
      node.chi_prime = tmp_dchi[i];
      node.alpha = TOA.adherence.alpha[i];
      eval_Omega_xyz(node);
      eval_g_xyz(node);
      if (this->m_ablation_flags.compute_Omega_xyz_prime)
      {
        eval_Omega_prime_xyz(node);
      }
      if (i > 0)
      {
        auto &cell = this->Cells[i - 1];
        cell.s_0 = TOA.reference.abscissa[i - 1];
        cell.s_1 = TOA.reference.abscissa[i];
        cell.L = cell.s_1 - cell.s_0;
        cell.ID0 = i - 1;
        cell.ID1 = i;
        cell.m_type = SegmentType::UNKNOWN;
      }
    }
  }

  static void eval_Omega_xyz(NodeStruct3D &node)
  {
    const real cos_mu = std::cos(node.mu);
    const real sin_mu = std::sin(node.mu);
    const real cos_phi = std::cos(node.phi);
    const real sin_phi = std::sin(node.phi);
    node.Omega_x = node.phi_prime - sin_mu * node.theta_prime;
    node.Omega_y = cos_mu * node.mu_prime + cos_mu * sin_phi * node.theta_prime;
    node.Omega_z = -sin_phi * node.mu_prime + cos_mu * cos_phi * node.theta_prime;
  }

  static void eval_Omega_prime_xyz(NodeStruct3D &node)
  {
    const real cos_mu = std::cos(node.mu);
    const real sin_mu = std::sin(node.mu);
    const real cos_phi = std::cos(node.phi);
    const real sin_phi = std::sin(node.phi);
    node.Omega_x_prime = node.phi_double_prime - cos_mu * node.mu_prime * node.theta_prime - sin_mu * node.theta_double_prime;
    const real t1 = cos_mu * node.theta_double_prime - node.mu_prime * (sin_mu * node.theta_prime + node.phi_prime);
    const real t2 = cos_mu * node.theta_prime * node.phi_prime + node.mu_double_prime;
    node.Omega_y_prime = t1 * sin_phi + cos_phi * t2;
    node.Omega_z_prime = t1 * cos_phi - sin_phi * t2;
  }

  static void eval_g_xyz(NodeStruct3D &node)
  {
    const real cos_mu = std::cos(node.mu);
    const real sin_mu = std::sin(node.mu);
    const real cos_phi = std::cos(node.phi);
    const real sin_phi = std::sin(node.phi);
    const real cos_chi = std::cos(node.chi);
    const real sin_chi = std::sin(node.chi);
    node.g_x = GRAVITY * (-sin_mu * cos_chi + cos_mu * sin_phi * sin_chi);
    node.g_y = GRAVITY * (sin_mu * sin_chi + cos_mu * sin_phi * cos_chi);
    node.g_z = GRAVITY * cos_mu * cos_phi;
  }

  // --------------------------------------------------------------------------------------------

  void compute_Vmax()
  {
    constexpr real k_small = 1e-6;
    // Deliberate deviation from FBGA_3D: FBGA_INDY/FBGA_MOTO both hardcoded this bracket to a
    // literal 130, so set_max_velocity() had no effect here at all (only on the initial-velocity
    // fallback in compute()). Confirmed and approved -- see FBGA3D_INTEGRATION_PLAN.md.
    const real v_top_speed = std::min(VMAX, this->max_velocity);
    real vmax = QUIET_NAN;
    for (auto &node : this->Nodes)
    {
      if (std::abs(node.Omega_z) < k_small)
      {
        vmax = v_top_speed;
      }
      else
      {
        constexpr real v_a = 0.0;
        const real v_b = v_top_speed;
        auto F2solve = [this, &node](const real v) -> real { return this->max_lateral_performance_func(node, v); };
        const bool ok = this->BD.solve(F2solve, v_a, v_b, vmax);
        vmax = ok ? vmax : v_top_speed;
      }
      node.V_max = vmax;
    }
  }

  [[nodiscard]] real max_lateral_performance_func(NodeStruct3D const &node, real V) const
  {
    const real V_dot = this->eval_V_dot_Vatildex(node, V, this->gggv.a_x_neutral(V));
    return std::abs(this->eval_a_tilde_y(node, V)) -
           (this->eval_a_tilde_y_lim(node, V, V_dot) + this->m_lat_tol - this->m_lat_tol_vmax) * this->m_lateral_shrink_factor;
  }

  [[nodiscard]] real get_max_longitudinal_performance(NodeStruct3D const &node, real V)
  {
    const NodeVelocityState st = this->compute_velocity_state(node, V);
    auto F2solve = [this, &node, &st, V](const real V_dot) -> real {
      return this->max_longitudinal_performance_func(node, st, V, V_dot);
    };

    const real V_dot_0 = eval_V_dot_Vatildex(node, st, this->gggv.a_x_neutral(V));
    // Self-consistent load-transfer estimate: both vehicles now use the converted V_dot_0 here
    // (FBGA_3D's FBGA_INDY originally reused the raw neutral value instead, which is only
    // equivalent on flat/non-rotating road segments -- see FBGA3D_INTEGRATION_PLAN.md).
    const real a_tilde_x_top = this->eval_a_tilde_x_max(node, st, V, V_dot_0);
    real V_dot_max = eval_V_dot_Vatildex(node, st, a_tilde_x_top);

    if ((std::abs(V_dot_max - V_dot_0) <= this->solver_p.tolerance) ||
        (std::abs(F2solve(V_dot_max)) <= this->solver_p.tolerance))
    {
      return V_dot_max;
    }
    if (std::abs(F2solve(V_dot_0)) <= this->solver_p.tolerance)
    {
      return V_dot_0;
    }
    while (F2solve(V_dot_max) < -this->solver_p.tolerance)
    {
      constexpr real scale_gain = 1.5;
      V_dot_max = (V_dot_max - V_dot_0) * scale_gain + V_dot_0;
    }
    real V_dot_sol = V_dot_0;
    const bool ok = this->BD.solve(F2solve, V_dot_0, V_dot_max, V_dot_sol);
    if (!ok)
    {
      std::cout << "Fbga3dSolver::get_max_longitudinal_performance: no solution for V = " << V << "\n";
    }
    return ok ? V_dot_sol : V_dot_0;
  }

  [[nodiscard]] real get_min_longitudinal_performance(NodeStruct3D const &node, real V)
  {
    const NodeVelocityState st = this->compute_velocity_state(node, V);
    auto F2solve = [this, &node, &st, V](const real V_dot) -> real {
      return this->min_longitudinal_performance_func(node, st, V, V_dot);
    };
    const real V_dot_0 = eval_V_dot_Vatildex(node, st, 0.0);
    const real a_tilde_x_bottom = this->eval_a_tilde_x_min(node, st, V, 0.0);
    real V_dot_min = eval_V_dot_Vatildex(node, st, a_tilde_x_bottom);
    if (std::abs(V_dot_min - V_dot_0) <= this->solver_p.tolerance)
    {
      return V_dot_min;
    }
    if (std::abs(F2solve(V_dot_min)) <= this->solver_p.tolerance)
    {
      return V_dot_min;
    }
    if (std::abs(F2solve(V_dot_0)) <= this->solver_p.tolerance)
    {
      return V_dot_0;
    }
    while (F2solve(V_dot_min) > +this->solver_p.tolerance)
    {
      V_dot_min = (V_dot_min - V_dot_0) * 1.5 + V_dot_0;
    }
    real V_dot_sol = V_dot_0;
    const bool ok = this->BD.solve(F2solve, V_dot_min, V_dot_0, V_dot_sol);
    if (!ok)
    {
      std::cout << "Fbga3dSolver::get_min_longitudinal_performance: no solution for V = " << V << "\n";
    }
    return ok ? V_dot_sol : V_dot_0;
  }

  [[nodiscard]] real min_longitudinal_performance_func(
    NodeStruct3D const &node, NodeVelocityState const &st, real V, real V_dot
  ) const
  {
    Context ctx;
    this->compute_context(node, st, V, V_dot, ctx);
    return ctx.a_tilde_x - ctx.a_tilde_x_min - this->solver_p.tolerance;
  }

  [[nodiscard]] real max_longitudinal_performance_func(
    NodeStruct3D const &node, NodeVelocityState const &st, real V, real V_dot
  ) const
  {
    Context ctx;
    this->compute_context(node, st, V, V_dot, ctx);
    return ctx.a_tilde_x - ctx.a_tilde_x_max + this->solver_p.tolerance;
  }

  [[nodiscard]] real get_min_longitudinal_performance_YF(NodeStruct3D const &node, real V)
  {
    const NodeVelocityState st = this->compute_velocity_state(node, V);
    auto F2solve = [this, &node, &st, V](const real V_dot) -> real {
      return this->min_longitudinal_performance_func_YF(node, st, V, V_dot);
    };
    const real V_dot_0 = eval_V_dot_Vatildex(node, st, 0.0);
    const real a_tilde_x_bottom = this->eval_a_tilde_x_min_YF(node, st, V, 0.0);
    real V_dot_min = eval_V_dot_Vatildex(node, st, a_tilde_x_bottom);
    if (std::abs(V_dot_min - V_dot_0) <= this->solver_p.tolerance)
    {
      return V_dot_min;
    }
    if (std::abs(F2solve(V_dot_min)) <= this->solver_p.tolerance)
    {
      return V_dot_min;
    }
    if (std::abs(F2solve(V_dot_0)) <= this->solver_p.tolerance)
    {
      return V_dot_0;
    }
    while (F2solve(V_dot_min) > +this->solver_p.tolerance)
    {
      V_dot_min = (V_dot_min - V_dot_0) * 1.5 + V_dot_0;
    }
    real V_dot_sol = V_dot_0;
    const bool ok = this->BD.solve(F2solve, V_dot_min, V_dot_0, V_dot_sol);
    if (!ok)
    {
      std::cout << "Fbga3dSolver::get_min_longitudinal_performance_YF: no solution for V = " << V << "\n";
    }
    return ok ? V_dot_sol : V_dot_0;
  }

  [[nodiscard]] real min_longitudinal_performance_func_YF(
    NodeStruct3D const &node, NodeVelocityState const &st, real V, real V_dot
  ) const
  {
    return eval_a_tilde_x(node, st, V_dot) - this->eval_a_tilde_x_min_YF(node, st, V, V_dot) - this->solver_p.tolerance;
  }

  // --------------------------------------------------------------------------------------------

  void BY()
  {
    bool ok_by = false;
    real v_0 = std::min(this->v_I, this->Nodes[0].V_max);
    this->Nodes[0].V_max = v_0;
    for (integer ith = 0; ith < static_cast<integer>(this->Cells.size()); ith++)
    {
      auto &cell = this->Cells[ith];
      auto &node_0 = this->Nodes[cell.ID0];
      auto &node_1 = this->Nodes[cell.ID1];
      cell.V_0 = v_0;

      if (v_0 <= this->m_yellow_flag_data.v_des_max)
      {
        this->yellow_index = ith;
        ok_by = true;
        break;
      }

      const real V_dot_0 = this->eval_V_dot_Vatildex(node_0, v_0, 0.0);
      real V_dot_min = this->get_min_longitudinal_performance_YF(node_0, v_0);
      const real a_hat_x_v_ref =
        ((this->m_yellow_flag_data.v_des_max * this->m_yellow_flag_data.v_des_max) - (v_0 * v_0)) / (2.0 * cell.L);
      const real V_dot_v_ref = this->eval_V_dot_Vahatx(node_0, v_0, a_hat_x_v_ref);

      V_dot_min = std::min(V_dot_min, V_dot_0);
      V_dot_min = std::max(V_dot_min, V_dot_v_ref);

      const real distance_V_dot_min_1 =
        this->signed_distance_YF(node_1, this->eval_V_next(node_0, cell, V_dot_min), V_dot_min);

      real VDOTSOL = 0.0;
      if (distance_V_dot_min_1 <= 0)
      {
        VDOTSOL = V_dot_min;
      }
      else
      {
        this->BD.solve(
          [this, &node_0, &node_1, &cell](const real VDOT) -> real {
            return this->signed_distance_YF(node_1, this->eval_V_next(node_0, cell, VDOT), VDOT) + this->solver_p.tolerance;
          },
          V_dot_min,
          V_dot_0,
          VDOTSOL
        );
      }
      if (!std::isnan(VDOTSOL))
      {
        cell.V_dot = VDOTSOL;
        cell.V_1 = this->eval_V_next(node_0, cell, VDOTSOL);
        v_0 = cell.V_1;
        cell.m_type = SegmentType::YELLOWFLAG;
        node_1.V_max = std::min(node_1.V_max, cell.V_1);
      }
      else
      {
        cell.m_type = SegmentType::YELLOWFLAG_NAN;
        this->dump_seg_id.push_back(cell.ID0);
        node_1.V_max = std::min(node_1.V_max, v_0);
        v_0 = node_1.V_max;
        cell.V_1 = v_0;
        this->yellow_index = ith;
        ok_by = false;
        break;
      }
    }
    for (integer ith = this->yellow_index; ith < static_cast<integer>(this->Cells.size()); ith++)
    {
      auto &cell = this->Cells[ith];
      auto &node_1 = this->Nodes[cell.ID1];
      if (node_1.V_max <= this->m_yellow_flag_data.v_des_max)
      {
        ok_by = true;
      }
      if (ok_by)
      {
        node_1.V_max = std::min(node_1.V_max, this->m_yellow_flag_data.v_des_max);
      }
    }
  }

  void FW()
  {
    real v_I = this->Cells[this->yellow_index].V_0;
    if (this->yellow_index == 0)
    {
      v_I = this->v_I;
    }
    real v_0 = std::min(v_I, this->Nodes[this->yellow_index].V_max);
    this->Nodes[this->yellow_index].V_max = v_0;
    for (integer ith = this->yellow_index; ith < static_cast<integer>(this->Cells.size()); ith++)
    {
      auto &cell = this->Cells[ith];
      auto &node_0 = this->Nodes[cell.ID0];
      auto &node_1 = this->Nodes[cell.ID1];
      cell.V_0 = v_0;

      const real a_hat_x_Vmax = ((node_1.V_max * node_1.V_max) - (v_0 * v_0)) / (2.0 * cell.L);
      const real V_dot_x_Vmax = this->eval_V_dot_Vahatx(node_0, v_0, a_hat_x_Vmax);
      const real V_dot_0 = this->eval_V_dot_Vatildex(node_0, v_0, 0.0);

      real V_dot_max = this->get_max_longitudinal_performance(node_0, v_0);
      V_dot_max = std::min(V_dot_max, V_dot_x_Vmax);

      const real distance_V_dot_max_1 = this->signed_distance(node_1, this->eval_V_next(node_0, cell, V_dot_max), V_dot_max);

      real VDOTSOL = 0.0;
      if (distance_V_dot_max_1 <= 0)
      {
        VDOTSOL = V_dot_max;
      }
      else
      {
        this->BD.solve(
          [this, &node_0, &node_1, &cell](const real VDOT) -> real {
            return this->signed_distance(node_1, this->eval_V_next(node_0, cell, VDOT), VDOT) + this->solver_p.tolerance;
          },
          V_dot_0,
          V_dot_max,
          VDOTSOL
        );
      }
      if (!std::isnan(VDOTSOL))
      {
        cell.V_dot = VDOTSOL;
        cell.V_1 = this->eval_V_next(node_0, cell, VDOTSOL);
        v_0 = cell.V_1;
        cell.m_type = SegmentType::FORWARD;
        node_1.V_max = std::min(node_1.V_max, cell.V_1);
      }
      else
      {
        cell.m_type = SegmentType::FORWARD_NAN;
        this->dump_seg_id.push_back(cell.ID0);
        node_1.V_max = std::min(node_1.V_max, v_0);
        v_0 = node_1.V_max;
        cell.V_1 = v_0;
      }
    }
  }

  void BW()
  {
    real v_1 = this->Nodes.back().V_max;
    if (!std::isnan(this->Cells.back().V_dot))
    {
      v_1 = this->Cells.back().V_1;
    }
    for (auto i = static_cast<integer>(this->Cells.size() - 1); i >= 0; i--)
    {
      auto &cell = this->Cells[i];
      auto &node_0 = this->Nodes[cell.ID0];
      auto &node_1 = this->Nodes[cell.ID1];
      cell.V_1 = v_1;

      const real V_dot_max = this->get_max_longitudinal_performance(node_1, v_1);
      const real V_dot_min = this->get_min_longitudinal_performance(node_1, v_1);

      const real v0_reach_max = std::min(node_0.V_max, this->eval_V_prev(node_1, cell, V_dot_min));
      const real v0_reach_min = std::max(0.0, this->eval_V_prev(node_1, cell, V_dot_max));

      const bool is_v0_reachable = (cell.V_0 >= v0_reach_min && cell.V_0 <= v0_reach_max);
      const real a_mean = (cell.V_1 * cell.V_1 - cell.V_0 * cell.V_0) / (2.0 * cell.L);
      const real V_dot_mean = this->eval_V_dot_Vahatx(node_0, cell.V_0, a_mean);
      const bool is_V_dot_mean_candidate = ((V_dot_mean >= V_dot_min) && (V_dot_mean <= V_dot_max));
      const real distance_0_mean = this->signed_distance(node_0, this->eval_V_prev(node_0, cell, V_dot_mean), V_dot_mean);
      const bool is_valid_forward =
        (cell.m_type == SegmentType::FORWARD &&
         (std::abs(this->eval_V_next(node_0, cell, cell.V_dot) - v_1) <= this->solver_p.tolerance));
      const bool is_valid_yellow_flag =
        (cell.m_type == SegmentType::YELLOWFLAG &&
         (std::abs(this->eval_V_next(node_0, cell, cell.V_dot) - v_1) <= this->solver_p.tolerance));
      if (is_valid_forward || is_valid_yellow_flag)
      {
        v_1 = cell.V_0;
        continue;
      }
      if (is_V_dot_mean_candidate && is_v0_reachable && (distance_0_mean <= (2 * this->solver_p.tolerance)))
      {
        cell.V_dot = V_dot_mean;
        cell.m_type = SegmentType::TRANSITION;
        v_1 = cell.V_0;
        continue;
      }
      else
      {
        const real distance_amin = this->signed_distance(node_0, this->eval_V_prev(node_0, cell, V_dot_min), V_dot_min);
        real VDOTSOL = 0.0;
        if (distance_amin <= 0)
        {
          VDOTSOL = V_dot_min;
        }
        else
        {
          const real V_dot_solver = this->eval_V_dot_Vatildex(node_1, cell.V_1, 0.0);
          this->BD.solve(
            [this, &node_0, &cell](const real VDOT) -> real {
              return this->signed_distance(node_0, this->eval_V_prev(node_0, cell, VDOT), VDOT) + this->solver_p.tolerance;
            },
            V_dot_min,
            V_dot_solver,
            VDOTSOL
          );
        }
        if (!std::isnan(VDOTSOL))
        {
          cell.V_dot = VDOTSOL;
          cell.V_0 = this->eval_V_prev(node_0, cell, VDOTSOL);
          cell.m_type = SegmentType::BACKWARD;
          v_1 = cell.V_0;
        }
        else
        {
          cell.V_dot = QUIET_NAN;
          cell.m_type = SegmentType::BACKWARD_NAN;
          this->dump_seg_id.push_back(cell.ID0);
        }
      }
    }
  }

  [[nodiscard]] real compute_total_time() const
  {
    real T = 0;
    for (const auto &cell : this->Cells)
    {
      T += static_cast<real>(2) * cell.L / (cell.V_0 + cell.V_1);
    }
    return T;
  }

  // --------------------------------------------------------------------------------------------

  void compute_context(NodeStruct3D const &node, real V, real V_dot, Context &ctx) const
  {
    this->compute_context(node, this->compute_velocity_state(node, V), V, V_dot, ctx);
  }

  void compute_context(NodeStruct3D const &node, NodeVelocityState const &st, real V, real V_dot, Context &ctx) const
  {
    const real s_ddotB = (st.chi_dot * st.cos_chi * V - st.cos_chi * V_dot) * st.inv_den_common;
    const real s_ddot = (!this->m_ablation_flags.compute_s_ddot) ? 0.0 : st.s_ddotA + s_ddotB;
    const real w_dot = node.Omega_x * st.s_dot * st.n_dot + node.Omega_x * s_ddot * node.n +
                        node.Omega_x_prime * st.s_dot * st.s_dot * node.n;

    ctx.a_tilde_x = (st.omega_hat_y * st.w) + V_dot + node.g_x;
    ctx.a_tilde_y = st.a_tilde_y;
    ctx.a_tilde_z = w_dot - (st.omega_hat_y * V) + node.g_z;

    ctx.a_tilde_y_lim = (node.alpha * this->gggv.a_y_lim(V, ctx.a_tilde_z) - this->m_lat_tol);
    ctx.a_tilde_y_clip = clip(ctx.a_tilde_y, -ctx.a_tilde_y_lim, ctx.a_tilde_y_lim);

    ctx.a_tilde_x_max_gg = this->gggv.a_x_push(ctx.a_tilde_y_clip, V, ctx.a_tilde_z, node.alpha);
    ctx.a_tilde_x_min_gg = this->gggv.a_x_pull(ctx.a_tilde_y_clip, V, ctx.a_tilde_z, node.alpha);
    ctx.a_tilde_x_eng = this->gggv.a_x_eng(V);

    ctx.a_tilde_x_max = std::min(ctx.a_tilde_x_max_gg, ctx.a_tilde_x_eng);
    ctx.a_tilde_x_min = ctx.a_tilde_x_min_gg;

    ctx.at_lateral_limit = (std::abs(ctx.a_tilde_y_clip) >= ctx.a_tilde_y_lim);
    if (ctx.at_lateral_limit)
    {
      ctx.a_tilde_x_max = std::min(0.0, ctx.a_tilde_x_eng);
      ctx.a_tilde_x_min = std::max(0.0, ctx.a_tilde_x_min);
    }
    if ((std::abs(ctx.a_tilde_y) - ctx.a_tilde_y_lim) > this->solver_p.tolerance)
    {
      ctx.a_tilde_x_max = std::min(0.0, ctx.a_tilde_x_eng);
      ctx.a_tilde_x_min = std::max(0.0, ctx.a_tilde_x_min);
    }
  }

  [[nodiscard]] real signed_distance(NodeStruct3D const &node, real V, real V_dot) const
  {
    Context ctx;
    this->compute_context(node, V, V_dot, ctx);
    return fb::utils::signed_distance(
      ctx.a_tilde_x, ctx.a_tilde_x_min, ctx.a_tilde_x_max, ctx.a_tilde_y, -ctx.a_tilde_y_lim, ctx.a_tilde_y_lim
    );
  }

  [[nodiscard]] real signed_distance_YF(NodeStruct3D const &node, real V, real V_dot) const
  {
    const real a_tilde_y = this->eval_a_tilde_y(node, V);
    const real a_tilde_lim_y = (this->eval_a_tilde_y_lim(node, V, V_dot) - this->m_lat_tol);
    const real a_tilde_x_max = this->eval_a_tilde_x_max(node, V, V_dot);
    const real a_tilde_x_min = this->eval_a_tilde_x_min_YF(node, V, V_dot);
    const real a_tilde_x = this->eval_a_tilde_x(node, V, V_dot);
    return fb::utils::signed_distance(a_tilde_x, a_tilde_x_min, a_tilde_x_max, a_tilde_y, -a_tilde_lim_y, a_tilde_lim_y);
  }

  [[nodiscard]] static real eval_V_next(NodeStruct3D const &node, CellStruct3D const &cell, real VDOT)
  {
    const real a_hat_x = eval_a_hat_x(node, cell.V_0, VDOT);
    return std::sqrt((static_cast<real>(2) * cell.L * a_hat_x) + std::pow(cell.V_0, 2));
  }

  [[nodiscard]] static real eval_V_prev(NodeStruct3D const &node, CellStruct3D const &cell, real VDOT)
  {
    const real a_hat_x = eval_a_hat_x(node, cell.V_1, VDOT);
    return std::sqrt(-static_cast<real>(2) * cell.L * a_hat_x + std::pow(cell.V_1, 2));
  }

  [[nodiscard]] static real eval_V(NodeStruct3D const &node, CellStruct3D const &cell, real s)
  {
    const real a_hat_x = eval_a_hat_x(node, cell.V_0, cell.V_dot);
    return std::sqrt((static_cast<real>(2) * s * a_hat_x) + std::pow(cell.V_0, 2));
  }

  [[nodiscard]] static real eval_a_tilde_x(NodeStruct3D const &node, real V, real V_dot)
  {
    return eval_a_hat_x(node, V, V_dot) + node.g_x;
  }

  [[nodiscard]] static real eval_a_tilde_x(NodeStruct3D const &node, NodeVelocityState const &st, real V_dot)
  {
    return eval_a_hat_x(st, V_dot) + node.g_x;
  }

  [[nodiscard]] static real eval_a_hat_x(NodeStruct3D const &node, real V, real V_dot)
  {
    const real cos_chi = std::cos(node.chi);
    const real sin_chi = std::sin(node.chi);
    const real s_dot = V * cos_chi / (1 - node.n * node.Omega_z);
    const real w = node.Omega_x * s_dot * node.n;
    const real omega_hat_y = (cos_chi * node.Omega_y - sin_chi * node.Omega_x) * s_dot;
    return (omega_hat_y * w) + V_dot;
  }

  [[nodiscard]] static real eval_a_hat_x(NodeVelocityState const &st, real V_dot)
  {
    return (st.omega_hat_y * st.w) + V_dot;
  }

  [[nodiscard]] static real eval_V_dot_Vatildex(NodeStruct3D const &node, real V, real atildex)
  {
    const real cos_chi = std::cos(node.chi);
    const real sin_chi = std::sin(node.chi);
    const real s_dot = V * cos_chi / (1 - node.n * node.Omega_z);
    const real w = node.Omega_x * s_dot * node.n;
    const real omega_hat_y = (cos_chi * node.Omega_y - sin_chi * node.Omega_x) * s_dot;
    return atildex - node.g_x - (omega_hat_y * w);
  }

  [[nodiscard]] static real eval_V_dot_Vatildex(NodeStruct3D const &node, NodeVelocityState const &st, real atildex)
  {
    return atildex - node.g_x - (st.omega_hat_y * st.w);
  }

  [[nodiscard]] static real eval_V_dot_Vahatx(NodeStruct3D const &node, real V, real ahatx)
  {
    const real cos_chi = std::cos(node.chi);
    const real sin_chi = std::sin(node.chi);
    const real s_dot = V * cos_chi / (1 - node.n * node.Omega_z);
    const real w = node.Omega_x * s_dot * node.n;
    const real omega_hat_y = (cos_chi * node.Omega_y - sin_chi * node.Omega_x) * s_dot;
    return ahatx - (omega_hat_y * w);
  }

  [[nodiscard]] real eval_a_tilde_y(NodeStruct3D const &node, real V) const
  {
    return this->eval_a_hat_y(node, V) + node.g_y;
  }

  [[nodiscard]] real eval_a_tilde_z(NodeStruct3D const &node, real V, real V_dot) const
  {
    return this->eval_a_hat_z(node, V, V_dot) + node.g_z;
  }

  [[nodiscard]] real eval_a_tilde_z(NodeStruct3D const &node, NodeVelocityState const &st, real V, real V_dot) const
  {
    return this->eval_a_hat_z(node, st, V, V_dot) + node.g_z;
  }

  [[nodiscard]] real eval_a_hat_y(NodeStruct3D const &node, real V) const
  {
    const real cos_chi = std::cos(node.chi);
    const real sin_chi = std::sin(node.chi);
    const real s_dot = V * cos_chi / (1 - node.n * node.Omega_z);
    const real w = node.Omega_x * s_dot * node.n;
    const real chi_dot = this->m_ablation_flags.compute_chi_dot ? node.chi_prime * s_dot : 0.0;
    const real omega_hat_x = (sin_chi * node.Omega_y + cos_chi * node.Omega_x) * s_dot;
    const real omega_hat_z = node.Omega_z * s_dot + chi_dot;
    return (omega_hat_z * V) - (omega_hat_x * w);
  }

  [[nodiscard]] real eval_a_hat_z(NodeStruct3D const &node, real V, real V_dot) const
  {
    const real cos_chi = std::cos(node.chi);
    const real sin_chi = std::sin(node.chi);
    const real s_dot = V * cos_chi / (1 - node.n * node.Omega_z);
    const real omega_hat_y = (cos_chi * node.Omega_y - sin_chi * node.Omega_x) * s_dot;
    const real n_dot = V * sin_chi;
    const real s_ddot = this->eval_s_ddot(node, V, V_dot);
    const real w_dot =
      node.Omega_x * s_dot * n_dot + node.Omega_x * s_ddot * node.n + node.Omega_x_prime * s_dot * s_dot * node.n;
    return w_dot - (omega_hat_y * V);
  }

  [[nodiscard]] real eval_a_hat_z(NodeStruct3D const &node, NodeVelocityState const &st, real V, real V_dot) const
  {
    const real s_ddotB = (st.chi_dot * st.cos_chi * V - st.cos_chi * V_dot) * st.inv_den_common;
    const real s_ddot = this->m_ablation_flags.compute_s_ddot ? st.s_ddotA + s_ddotB : 0.0;
    const real w_dot = node.Omega_x * st.s_dot * st.n_dot + node.Omega_x * s_ddot * node.n +
                        node.Omega_x_prime * st.s_dot * st.s_dot * node.n;
    return w_dot - (st.omega_hat_y * V);
  }

  [[nodiscard]] real eval_a_tilde_y_lim(NodeStruct3D const &node, real V, real V_dot) const
  {
    const real a_tilde_z = this->eval_a_tilde_z(node, V, V_dot);
    return node.alpha * this->gggv.a_y_lim(V, a_tilde_z);
  }

  [[nodiscard]] real eval_a_tilde_x_max(NodeStruct3D const &node, real V, real V_dot) const
  {
    const real a_tilde_y = this->eval_a_tilde_y(node, V);
    const real a_tilde_z = this->eval_a_tilde_z(node, V, V_dot);
    return this->gggv.a_x_push(a_tilde_y, V, a_tilde_z, node.alpha);
  }

  [[nodiscard]] real eval_a_tilde_x_max(NodeStruct3D const &node, NodeVelocityState const &st, real V, real V_dot) const
  {
    const real a_tilde_z = this->eval_a_tilde_z(node, st, V, V_dot);
    return this->gggv.a_x_push(st.a_tilde_y, V, a_tilde_z, node.alpha);
  }

  [[nodiscard]] real eval_a_tilde_x_min(NodeStruct3D const &node, real V, real V_dot) const
  {
    const real a_tilde_y = this->eval_a_tilde_y(node, V);
    const real a_tilde_z = this->eval_a_tilde_z(node, V, V_dot);
    return this->gggv.a_x_pull(a_tilde_y, V, a_tilde_z, node.alpha);
  }

  [[nodiscard]] real eval_a_tilde_x_min(NodeStruct3D const &node, NodeVelocityState const &st, real V, real V_dot) const
  {
    const real a_tilde_z = this->eval_a_tilde_z(node, st, V, V_dot);
    return this->gggv.a_x_pull(st.a_tilde_y, V, a_tilde_z, node.alpha);
  }

  [[nodiscard]] real eval_a_tilde_x_min_YF(NodeStruct3D const &node, real V, real V_dot) const
  {
    const real a_tilde_x_min = this->eval_a_tilde_x_min(node, V, V_dot);
    return std::max(a_tilde_x_min, this->m_yellow_flag_data.a_des_min);
  }

  [[nodiscard]] real eval_a_tilde_x_min_YF(NodeStruct3D const &node, NodeVelocityState const &st, real V, real V_dot) const
  {
    const real a_tilde_x_min = this->eval_a_tilde_x_min(node, st, V, V_dot);
    return std::max(a_tilde_x_min, this->m_yellow_flag_data.a_des_min);
  }

  [[nodiscard]] real eval_s_ddot(NodeStruct3D const &node, real V, real V_dot) const
  {
    if (!this->m_ablation_flags.compute_s_ddot)
    {
      return 0.0;
    }
    const real cos_chi = std::cos(node.chi);
    const real sin_chi = std::sin(node.chi);
    const real den_common = node.Omega_z * node.n - 1.0;
    const real inv_den_common = 1.0 / den_common;
    const real inv_den_common2 = inv_den_common * inv_den_common;
    const real s_dot = V * cos_chi * (-inv_den_common);
    const real chi_dot = this->m_ablation_flags.compute_chi_dot ? node.chi_prime * s_dot : 0.0;
    const real s_ddotA = -V * cos_chi * (node.Omega_z_prime * s_dot * node.n + node.Omega_z * V * sin_chi) * inv_den_common2;
    const real s_ddotB = (chi_dot * cos_chi * V - cos_chi * V_dot) * inv_den_common;
    return s_ddotA + s_ddotB;
  }
};

} // namespace fb::fbga3d
