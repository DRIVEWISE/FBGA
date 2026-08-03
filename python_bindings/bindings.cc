#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>

#include <utils/types.hh>
#include <fbga/fb2d.hh>

namespace py = pybind11;
using namespace fb::utils;
using namespace fb::fbga;

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

  // Add some useful constants
  m.attr("GRAVITY") = GRAVITY;
  m.attr("PI") = PI;
  m.attr("DEG2RAD") = DEG2RAD;
  m.attr("RAD2DEG") = RAD2DEG;
}
