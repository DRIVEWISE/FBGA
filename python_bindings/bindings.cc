#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>

#include "FBGA/types.hxx"
#include "FBGA/FWBW.hxx"

namespace py = pybind11;
using namespace GG;

PYBIND11_MODULE(fbga_py, m) {
  m.doc() = "Python bindings for the Forward-Backward with Generic Acceleration constraint (FBGA) library";

  // Bind the gg_range_max_min struct
  py::class_<gg_range_max_min>(m, "GGRangeMaxMin")
      .def(py::init<>())
      .def_readwrite("min", &gg_range_max_min::min, "Minimum function")
      .def_readwrite("max", &gg_range_max_min::max, "Maximum function");

  // Bind the FWBW class (essential methods only)
  py::class_<FWBW>(m, "FWBW")
      .def(py::init<const std::function<real(real, real)>&,
                     const std::function<real(real, real)>&,
                     const gg_range_max_min&>(),
           "Constructor with function pointers",
           py::arg("gg_Upper"), py::arg("gg_Lower"), py::arg("gg_range"))
      .def("compute", &FWBW::compute,
           "Compute the forward-backward algorithm",
           py::arg("SS"), py::arg("KK"), py::arg("v0"), py::arg("vfmax") = VMAX_SPEED)
      .def("compute_timing", &FWBW::compute_timing,
           "compute_timing the forward-backward algorithm",
           py::arg("SS"), py::arg("KK"), py::arg("v0"), py::arg("vfmax") = VMAX_SPEED)
      .def("evaluate",
           [](FWBW &self, const std::vector<real> &SS) {
             // FWBW::evaluate writes into AX/AY/V via operator[] without
             // resizing them itself, so the caller must pre-size them.
             std::vector<real> AX(SS.size()), AY(SS.size()), V(SS.size());
             self.evaluate(SS, AX, AY, V);
             return std::make_tuple(AX, AY, V);
           },
           "Evaluate acceleration and velocity at given positions; returns (AX, AY, V)",
           py::arg("SS"))
      .def("evalV", &FWBW::evalV,
           "Evaluate velocity at position s",
           py::arg("s"))
      .def("evalAx", &FWBW::evalAx,
           "Evaluate longitudinal acceleration at position s",
           py::arg("s"))
      .def("evalAy", &FWBW::evalAy,
           "Evaluate lateral acceleration at position s",
           py::arg("s"))
      .def("get_seg_idx", &FWBW::get_seg_idx,
           "Get segment index for position s",
           py::arg("s"))
      .def("evalVmax", &FWBW::evalVmax,
           "Evaluate maximum velocity at position s",
           py::arg("s"))
      .def("evalS", &FWBW::evalS,
           "Evaluate position at time t",
           py::arg("t"))
      .def("evalV_t", &FWBW::evalV_t,
           "Evaluate velocity at time t",
           py::arg("t"))
      .def("evalAx_t", &FWBW::evalAx_t,
           "Evaluate acc x at time t",
           py::arg("t"))
      .def("evalAy_t", &FWBW::evalAy_t,
           "Evaluate acc y at time t",
           py::arg("t"))
      .def("evalSegmentType", &FWBW::evalSegmentType,
           "Evaluate segment type at time t",
           py::arg("t"));

  // Add some useful constants
  m.attr("GRAVITY") = GRAVITY;
  m.attr("PI") = PI;
  m.attr("DEG2RAD") = DEG2RAD;
  m.attr("RAD2DEG") = RAD2DEG;
}