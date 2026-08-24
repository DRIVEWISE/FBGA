/*
(***********************************************************************************)
(*                                                                                 *)
(* The FBGA project                                                               *)
(*                                                                                 *)
(* Copyright (c) 2025, Mattia Piazza                                               *)
(*                                                                                 *)
(*    Mattia Piazza                                                                *)
(*    Department of Industrial Engineering                                         *)
(*    University of Trento                                                         *)
(*    e-mail: mattia.piazza@unitn.it                                               *)
(*                                                                                 *)
(***********************************************************************************)
*/

#pragma once

#include <utils/types.hh>
#include <utils/segment.hh>
#include <solvers/brent_dekker.hh>
#include <functional>

namespace fb::fbga
{

using fb::utils::real;
using fb::utils::integer;
using fb::utils::Segment;
using fb::solvers::BrentDekker;

using GgRangeMaxMin = struct GgRangeMaxMin
{
  std::function<real(real)> min = nullptr;
  std::function<real(real)> max = nullptr;
};

using SolverParams = struct SolverParams
{
  real tolerance        = STD_TOL;
  int max_iter          = STD_MAX_ITER;
  std::string verbosity = STD_VERBOSE;
};

using Solution2D = struct Solution2D
{
  std::vector<real> S;
  std::vector<real> AX;
  std::vector<real> AY;
  std::vector<real> V;
  std::vector<real> T;
};

#define VMAX_SPEED 130.0

class Fb2d
{
private:
  /* data */
  std::function<real(real, real)> gg_Upper = nullptr; // Upper bound function
  std::function<real(real, real)> gg_Lower = nullptr; // Lower bound function
  GgRangeMaxMin gg_range;         // Range of the curvature (struct with min and max functions)
  BrentDekker BD;                 // Solver
  SolverParams solver_p;          // Solver parameters
  std::vector<Segment> Segments; // Vector of segments
  std::vector<real> Vmax_vec;    // Vector of maximum reachable velocities
  std::vector<real> S_vec;       // Vector of abscissas
  std::vector<real> K_vec;       // Vector of curvatures
  real v_I{0.0};                 // Initial velocity
  real v_F{0.0};                 // Final velocity
  std::vector<int> dump_seg_id;  // Vector of segments with problems for debug
  real max_velocity{VMAX_SPEED}; // maximum velocity allowedq
protected:
  Fb2d();
  void set_max_velocity(real velocity) {
    max_velocity = velocity;
  }
public:
  // constructors
  Fb2d(
    const std::function<real(real, real)> &gg_Upper,
    const std::function<real(real, real)> &gg_Lower,
    const GgRangeMaxMin &gg_range
  );
  void setup_functions(
    const std::function<real(real, real)> &gg_Upper,
    const std::function<real(real, real)> &gg_Lower,
    const GgRangeMaxMin &gg_range
  );
  // main methods
  // core Forward-Backward method
  real compute(std::vector<real> const &SS, std::vector<real> const &KK, real v0, real vfmax = VMAX_SPEED);
  real compute_cyclic(std::vector<real> const &SS, std::vector<real> const &KK);
  real compute_timing(std::vector<real> const &SS, std::vector<real> const &KK, real v0, real vfmax = VMAX_SPEED);
  private:
  // compute Vmax vector
  void compute_Vmax();
  // Forward step
  void FW();
  // Backward step
  void BW();
  // compute time
  [[nodiscard]] real compute_time();
  public:
  // compute the distance with sign.
  [[nodiscard]] real signed_distance(real ax, real ay, real v) const;
  // check if a point is in the range
  [[nodiscard]] bool is_in_range(real ax, real ay, real v) const;
  // evaluation
  void evaluate(
    std::vector<real> const &SS, std::vector<real> &AX, std::vector<real> &AY, std::vector<real> &V, std::vector<real> &TT
  );
  void evaluate( Solution2D &sol );
  void evaluate_t(
    std::vector<real> const &TT, std::vector<real> &AX, std::vector<real> &AY, std::vector<real> &V, std::vector<real> &SS
  );
  void evaluate_t( Solution2D &sol );
  [[nodiscard]] real evalV(real s) const;
  [[nodiscard]] real evalT(real s) const;
  [[nodiscard]] real evalAx(real s) const;
  [[nodiscard]] real evalAy(real s) const;
  [[nodiscard]] integer get_seg_idx(real s) const;
  [[nodiscard]] real evalVmax(const real s) const;
  [[nodiscard]] real evalS(real t) const;
  [[nodiscard]] integer get_seg_idx_t(real t) const;

  [[nodiscard]] real evalV_t(const real t) const;
  [[nodiscard]] real evalAx_t(const real t) const;
  [[nodiscard]] real evalAy_t(const real t) const;

  [[nodiscard]] real evalSegmentType(const real t) const;

  // get dump
  [[nodiscard]] std::vector<int> get_dump() const { return this->dump_seg_id; }

  void check_segments();
};

} // namespace fb::fbga
