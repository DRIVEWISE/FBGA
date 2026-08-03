#pragma once

#include <fbga3d/fbga3d_solver.hh>
#include <fbga3d/gggv_indy.hh>

namespace fb::fbga3d
{

using FbgaIndy = Fbga3dSolver<GggvIndy>;

} // namespace fb::fbga3d

extern template class fb::fbga3d::Fbga3dSolver<fb::fbga3d::GggvIndy>;
