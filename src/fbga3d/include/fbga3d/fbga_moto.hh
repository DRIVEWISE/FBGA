#pragma once

#include <fbga3d/fbga3d_solver.hh>
#include <fbga3d/gggv_moto.hh>

namespace fb::fbga3d
{

// FBGA_3D's FBGA_MOTO used tighter lateral tolerances than FBGA_INDY's (which default to 0.0
// via the primary Fbga3dDefaults template) -- see FBGA3D_INTEGRATION_PLAN.md.
template <>
struct Fbga3dDefaults<GggvMoto>
{
  static constexpr real lat_tol = 1e-6;
  static constexpr real lat_tol_vmax = 1e-4;
};

using FbgaMoto = Fbga3dSolver<GggvMoto>;

} // namespace fb::fbga3d

extern template class fb::fbga3d::Fbga3dSolver<fb::fbga3d::GggvMoto>;
