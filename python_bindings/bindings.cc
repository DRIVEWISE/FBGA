#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>

#include <utils/types.hh>
#include <utils/gg_utils.hh>
#include <fbga/fb2d.hh>
#include <fbga3d/fbga_indy.hh>
#include <fbga3d/fbga_moto.hh>

namespace py = pybind11;
using namespace fb::utils;
using namespace fb::fbga;
using namespace fb::fbga3d;

namespace
{

// FbgaIndy and FbgaMoto are both Fbga3dSolver<GggvModel> instantiations with an identical public
// API -- bind the shared surface once instead of duplicating ~20 .def() calls per vehicle.
template <class GggvModel>
void bind_fbga3d_solver(py::module_ &m, char const *name)
{
  using Solver = fb::fbga3d::Fbga3dSolver<GggvModel>;
  py::class_<Solver>(m, name)
      .def(py::init<>())
      .def(
          "model",
          static_cast<GggvModel &(Solver::*)()>(&Solver::model),
          "The underlying GG-diagram model (mutable reference)",
          py::return_value_policy::reference_internal
      )
      .def("set_max_velocity", &Solver::set_max_velocity, py::arg("velocity"))
      .def(
          "compute", &Solver::compute, "Compute the velocity profile", py::arg("TOA"), py::arg("V0") = -1.0,
          py::arg("yellow_flag") = YellowFlagData()
      )
      .def("get_solution", &Solver::get_solution, "Get the computed solution")
      .def("get_dump", &Solver::get_dump, "Indices of cells that failed to converge during compute()")
      .def("get_cell_idx", &Solver::get_cell_idx, py::arg("s"))
      .def("eval_Vmax", &Solver::eval_Vmax, py::arg("s"))
      .def("eval_Omega_x", &Solver::eval_Omega_x, py::arg("s"))
      .def("eval_Omega_y", &Solver::eval_Omega_y, py::arg("s"))
      .def("eval_Omega_z", &Solver::eval_Omega_z, py::arg("s"))
      // eval_V is overloaded (this public accessor vs. a private static helper with a different
      // signature) -- disambiguate with an explicit cast, same as model() above.
      .def("eval_V", static_cast<real (Solver::*)(real) const>(&Solver::eval_V), py::arg("s"))
      .def("eval_V_dot", &Solver::eval_V_dot, py::arg("s"))
      .def("eval_A_hat_x", &Solver::eval_A_hat_x, py::arg("s"))
      .def("eval_A_hat_y", &Solver::eval_A_hat_y, py::arg("s"))
      .def("eval_A_hat_z", &Solver::eval_A_hat_z, py::arg("s"))
      .def("eval_A_tilde_x", &Solver::eval_A_tilde_x, py::arg("s"))
      .def("eval_A_tilde_y", &Solver::eval_A_tilde_y, py::arg("s"))
      .def("eval_A_tilde_z", &Solver::eval_A_tilde_z, py::arg("s"))
      .def("eval_g_x", &Solver::eval_g_x, py::arg("s"))
      .def("eval_g_y", &Solver::eval_g_y, py::arg("s"))
      .def("eval_g_z", &Solver::eval_g_z, py::arg("s"))
      .def("eval_segment_type", &Solver::eval_segment_type, py::arg("s"))
      .def("eval_alpha", &Solver::eval_alpha, py::arg("s"));
}

} // namespace

PYBIND11_MODULE(fbga_py, m) {
  m.doc() = "Python bindings for the Forward-Backward with Generic Acceleration constraint (FBGA) library";

  // Bind the GgRangeMaxMin struct
  py::class_<GgRangeMaxMin>(m, "GgRangeMaxMin")
      .def(py::init<>())
      .def_readwrite("min", &GgRangeMaxMin::min, "Minimum function")
      .def_readwrite("max", &GgRangeMaxMin::max, "Maximum function");

  // Bind the Fb2d class (essential methods only)
  py::class_<Fb2d>(m, "Fb2d")
      .def(py::init<const std::function<real(real, real)>&,
                     const std::function<real(real, real)>&,
                     const GgRangeMaxMin&>(),
           "Constructor with function pointers",
           py::arg("gg_Upper"), py::arg("gg_Lower"), py::arg("gg_range"))
      .def("compute", &Fb2d::compute,
           "Compute the forward-backward algorithm",
           py::arg("SS"), py::arg("KK"), py::arg("v0"), py::arg("vfmax") = VMAX_SPEED)
      .def("compute_timing", &Fb2d::compute_timing,
           "compute_timing the forward-backward algorithm",
           py::arg("SS"), py::arg("KK"), py::arg("v0"), py::arg("vfmax") = VMAX_SPEED)
      .def("evaluate",
           [](Fb2d &self, const std::vector<real> &SS) {
             // Fb2d::evaluate writes into AX/AY/V via operator[] without
             // resizing them itself, so the caller must pre-size them.
             std::vector<real> AX(SS.size()), AY(SS.size()), V(SS.size());
             self.evaluate(SS, AX, AY, V);
             return std::make_tuple(AX, AY, V);
           },
           "Evaluate acceleration and velocity at given positions; returns (AX, AY, V)",
           py::arg("SS"))
      .def("evalV", &Fb2d::evalV,
           "Evaluate velocity at position s",
           py::arg("s"))
      .def("evalAx", &Fb2d::evalAx,
           "Evaluate longitudinal acceleration at position s",
           py::arg("s"))
      .def("evalAy", &Fb2d::evalAy,
           "Evaluate lateral acceleration at position s",
           py::arg("s"))
      .def("get_seg_idx", &Fb2d::get_seg_idx,
           "Get segment index for position s",
           py::arg("s"))
      .def("evalVmax", &Fb2d::evalVmax,
           "Evaluate maximum velocity at position s",
           py::arg("s"))
      .def("evalS", &Fb2d::evalS,
           "Evaluate position at time t",
           py::arg("t"))
      .def("evalV_t", &Fb2d::evalV_t,
           "Evaluate velocity at time t",
           py::arg("t"))
      .def("evalAx_t", &Fb2d::evalAx_t,
           "Evaluate acc x at time t",
           py::arg("t"))
      .def("evalAy_t", &Fb2d::evalAy_t,
           "Evaluate acc y at time t",
           py::arg("t"))
      .def("evalSegmentType", &Fb2d::evalSegmentType,
           "Evaluate segment type at time t",
           py::arg("t"));

  // -----------------------------------------------------------------------------------------
  // fbga3d: 3D forward-backward solver, generic over the vehicle's GG-diagram model.

  py::enum_<SegmentType>(m, "SegmentType")
      .value("UNKNOWN", SegmentType::UNKNOWN)
      .value("FORWARD", SegmentType::FORWARD)
      .value("BACKWARD", SegmentType::BACKWARD)
      .value("TRANSITION", SegmentType::TRANSITION)
      .value("FORWARD_NAN", SegmentType::FORWARD_NAN)
      .value("BACKWARD_NAN", SegmentType::BACKWARD_NAN)
      .value("TRANSITION_NAN", SegmentType::TRANSITION_NAN)
      .value("YELLOWFLAG", SegmentType::YELLOWFLAG)
      .value("YELLOWFLAG_NAN", SegmentType::YELLOWFLAG_NAN);

  py::class_<ScalingGggvFactors>(m, "ScalingGggvFactors")
      .def(py::init<>())
      .def(py::init<real, real, real, real, real>(),
           py::arg("ax_max_scale"), py::arg("ax_min_scale"), py::arg("ay_scale"),
           py::arg("gg_exponent_ax_pos"), py::arg("gg_exponent_ax_neg"))
      .def_readwrite("ax_max_scale", &ScalingGggvFactors::ax_max_scale)
      .def_readwrite("ax_min_scale", &ScalingGggvFactors::ax_min_scale)
      .def_readwrite("ay_scale", &ScalingGggvFactors::ay_scale)
      .def_readwrite("gg_exponent_ax_pos", &ScalingGggvFactors::gg_exponent_ax_pos)
      .def_readwrite("gg_exponent_ax_neg", &ScalingGggvFactors::gg_exponent_ax_neg);

  py::class_<SplineDataCollection>(m, "SplineDataCollection")
      .def(py::init<>())
      .def(py::init<
               std::vector<real> const &, std::vector<real> const &, std::vector<real> const &,
               std::vector<real> const &, std::vector<real> const &, ScalingGggvFactors const &,
               std::vector<real> const &, std::vector<real> const &>(),
           py::arg("v_data"), py::arg("az_data"), py::arg("ay_max_data"), py::arg("ax_min_data"),
           py::arg("ax_max_data"), py::arg("scaling_factors"), py::arg("v_data_eng"), py::arg("ax_eng_data"))
      .def_readwrite("v_data", &SplineDataCollection::v_data)
      .def_readwrite("az_data", &SplineDataCollection::az_data)
      .def_readwrite("ay_max_data", &SplineDataCollection::ay_max_data)
      .def_readwrite("ax_min_data", &SplineDataCollection::ax_min_data)
      .def_readwrite("ax_max_data", &SplineDataCollection::ax_max_data)
      .def_readwrite("scaling_factors", &SplineDataCollection::scaling_factors)
      .def_readwrite("v_data_eng", &SplineDataCollection::v_data_eng)
      .def_readwrite("ax_eng_data", &SplineDataCollection::ax_eng_data);

  py::class_<EngineMaxGgv>(m, "EngineMaxGgv")
      .def(py::init<>())
      .def(py::init<std::vector<real> const &, std::vector<real> const &>(),
           py::arg("v_data_eng"), py::arg("ax_eng_data"))
      .def_readwrite("v_data_eng", &EngineMaxGgv::v_data_eng)
      .def_readwrite("ax_eng_data", &EngineMaxGgv::ax_eng_data);

  py::class_<YellowFlagData>(m, "YellowFlagData")
      .def(py::init<>())
      .def(py::init<real, real, bool>(), py::arg("v_des_max"), py::arg("a_des_min"), py::arg("is_yellow"))
      .def_readwrite("v_des_max", &YellowFlagData::v_des_max)
      .def_readwrite("a_des_min", &YellowFlagData::a_des_min)
      .def_readwrite("is_yellow", &YellowFlagData::is_yellow);

  py::class_<RoadAnglesAndDerivativesContainer>(m, "RoadAnglesAndDerivativesContainer")
      .def(py::init<>())
      .def(
          py::init(
              [](std::vector<real> const &mu, std::vector<real> const &phi, std::vector<real> const &theta,
                 std::vector<real> const &mu_prime, std::vector<real> const &phi_prime,
                 std::vector<real> const &theta_prime, std::vector<real> const &abscissa) {
                RoadAnglesAndDerivativesContainer c;
                c.mu = mu;
                c.phi = phi;
                c.theta = theta;
                c.mu_prime = mu_prime;
                c.phi_prime = phi_prime;
                c.theta_prime = theta_prime;
                c.abscissa = abscissa;
                return c;
              }
          ),
          py::arg("mu"), py::arg("phi"), py::arg("theta"), py::arg("mu_prime"), py::arg("phi_prime"),
          py::arg("theta_prime"), py::arg("abscissa")
      )
      .def_readwrite("mu", &RoadAnglesAndDerivativesContainer::mu)
      .def_readwrite("phi", &RoadAnglesAndDerivativesContainer::phi)
      .def_readwrite("theta", &RoadAnglesAndDerivativesContainer::theta)
      .def_readwrite("mu_prime", &RoadAnglesAndDerivativesContainer::mu_prime)
      .def_readwrite("phi_prime", &RoadAnglesAndDerivativesContainer::phi_prime)
      .def_readwrite("theta_prime", &RoadAnglesAndDerivativesContainer::theta_prime)
      .def_readwrite("abscissa", &RoadAnglesAndDerivativesContainer::abscissa);

  py::class_<TrajectoryOffsetContainer>(m, "TrajectoryOffsetContainer")
      .def(py::init<>())
      .def(
          py::init(
              [](std::vector<real> const &n, std::vector<real> const &chi) {
                TrajectoryOffsetContainer c;
                c.n = n;
                c.chi = chi;
                return c;
              }
          ),
          py::arg("n"), py::arg("chi")
      )
      .def_readwrite("n", &TrajectoryOffsetContainer::n)
      .def_readwrite("chi", &TrajectoryOffsetContainer::chi);

  py::class_<AdherenceContainer>(m, "AdherenceContainer")
      .def(py::init<>())
      .def(
          py::init(
              [](std::vector<real> const &alpha) {
                AdherenceContainer c;
                c.alpha = alpha;
                return c;
              }
          ),
          py::arg("alpha")
      )
      .def_readwrite("alpha", &AdherenceContainer::alpha);

  py::class_<TrajectoryOffsetAndAnglesContainer>(m, "TrajectoryOffsetAndAnglesContainer")
      .def(py::init<>())
      .def(
          py::init(
              [](TrajectoryOffsetContainer const &offset, RoadAnglesAndDerivativesContainer const &reference,
                 AdherenceContainer const &adherence) {
                TrajectoryOffsetAndAnglesContainer c;
                c.offset = offset;
                c.reference = reference;
                c.adherence = adherence;
                return c;
              }
          ),
          py::arg("offset"), py::arg("reference"), py::arg("adherence")
      )
      .def_readwrite("offset", &TrajectoryOffsetAndAnglesContainer::offset)
      .def_readwrite("reference", &TrajectoryOffsetAndAnglesContainer::reference)
      .def_readwrite("adherence", &TrajectoryOffsetAndAnglesContainer::adherence);

  py::class_<SolutionContainer>(m, "SolutionContainer")
      .def(py::init<>())
      .def_readwrite("V0", &SolutionContainer::V0)
      .def_readwrite("V1", &SolutionContainer::V1)
      .def_readwrite("V_dot", &SolutionContainer::V_dot);

  py::class_<GggvIndy>(m, "GggvIndy")
      .def(py::init<>(), "Default constructor; loads spline data from ./data/INDY/*.npy")
      .def(py::init<ScalingGggvFactors const &>(), py::arg("scaling_factors"))
      .def(py::init<SplineDataCollection const &>(), py::arg("spline_data"))
      .def("set_scaling_factors", &GggvIndy::set_scaling_factors, py::arg("scaling_factors"))
      .def("set_engine_max", &GggvIndy::set_engine_max, py::arg("engine_max_ggv"))
      .def("a_x_push", &GggvIndy::a_x_push, py::arg("ay_tilde"), py::arg("V"), py::arg("az_tilde"), py::arg("alpha") = 1.0)
      .def("a_x_pull", &GggvIndy::a_x_pull, py::arg("ay_tilde"), py::arg("V"), py::arg("az_tilde"), py::arg("alpha") = 1.0)
      .def("a_x_neutral", &GggvIndy::a_x_neutral, py::arg("V"))
      .def("rho", &GggvIndy::rho, py::arg("V"), py::arg("az_tilde"))
      .def("rho_max", &GggvIndy::rho_max, py::arg("V"), py::arg("az_tilde"))
      .def("rho_min", &GggvIndy::rho_min, py::arg("V"), py::arg("az_tilde"))
      .def("a_x_eng", &GggvIndy::a_x_eng, py::arg("V"))
      .def("a_x_max", &GggvIndy::a_x_max, py::arg("V"), py::arg("az_tilde"))
      .def("a_x_min", &GggvIndy::a_x_min, py::arg("V"), py::arg("az_tilde"))
      .def("a_y_lim", &GggvIndy::a_y_lim, py::arg("V"), py::arg("az_tilde"));

  py::class_<GggvMoto>(m, "GggvMoto")
      .def(py::init<>())
      .def("a_x_push", &GggvMoto::a_x_push, py::arg("ay_tilde"), py::arg("V"), py::arg("az_tilde"), py::arg("alpha") = 1.0)
      .def("a_x_pull", &GggvMoto::a_x_pull, py::arg("ay_tilde"), py::arg("V"), py::arg("az_tilde"), py::arg("alpha") = 1.0)
      .def("a_x_neutral", &GggvMoto::a_x_neutral, py::arg("V"))
      .def("a_x_eng", &GggvMoto::a_x_eng, py::arg("V"))
      .def("a_y_lim", &GggvMoto::a_y_lim, py::arg("V"), py::arg("az_tilde"))
      .def("a_x_aero", &GggvMoto::a_x_aero, py::arg("V"));

  // FbgaIndy/FbgaMoto: Fbga3dSolver<GggvIndy>/<GggvMoto> -- see fbga3d_solver.hh. Plotting
  // (eval_shell_plot*) and constraint-satisfaction-checking helpers from the original FBGA_3D
  // aren't ported to the C++ core yet, so they aren't bound here either.
  bind_fbga3d_solver<GggvIndy>(m, "FbgaIndy");
  bind_fbga3d_solver<GggvMoto>(m, "FbgaMoto");

  // Add some useful constants
  m.attr("GRAVITY") = GRAVITY;
  m.attr("PI") = PI;
  m.attr("DEG2RAD") = DEG2RAD;
  m.attr("RAD2DEG") = RAD2DEG;
}
