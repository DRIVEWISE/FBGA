#include <cmath>
#include <iostream>

#include <fbga/fb2d.hh>
#include <utils/types.hh>
#include <utils/gg_utils.hh>

using namespace fb::utils;
using namespace fb::fbga;

int main() {
  std::cout << "FBGA - Forward Backward Generic Acceleration constraints\n";
  std::cout << " > Running example: basic\n";
  // 1. Define track data
  std::vector<real> SS_vec = {0.0,50.0, 100.0, 150.0, 200.0};  // Arc length [m]
  std::vector<real> KK_vec = {0.001, 0.01, 0.005, -0.01, 0.0};    // Curvature [1/m]
  real v_initial = 20.0;  // Initial velocity [m/s]

  std::cout << " > Track data:\n";
  std::cout << " >   SS_vec: " << SS_vec.size() << " points\n";
  std::cout << " >   KK_vec: " << KK_vec.size() << " points\n";
  std::cout << " >   v_initial: " << v_initial << " m/s\n";
  for (size_t i = 0; i < SS_vec.size(); ++i) {
    std::cout << " >   SS[" << i << "]: " << SS_vec[i]
              << ", KK[" << i << "]: " << KK_vec[i] << "\n";
  }


  LinearInterpolatorSet spline_1D(
    SS_vec,
    { KK_vec },
    {"kappa"});

  integer numpt_eval = std::ceil(SS_vec.back());

  std::vector<real> SS_eval_vec(numpt_eval);
  std::vector<real> KK_eval_vec(numpt_eval);
  for (integer i = 0; i < std::ceil(SS_vec.back()); ++i) {
    SS_eval_vec[i] = static_cast<real>(i);
    KK_eval_vec[i] = spline_1D.eval("kappa", static_cast<real>(i));
  }


  // 2. Define vehicle parameters
  const real mu_x = 1.3;   // Longitudinal friction coefficient
  const real mu_y = 1.4;   // Lateral friction coefficient
  const real g = 9.81;     // Gravity [m/s²]

  // 3. Define constraint functions
  auto GG_shape1 = [mu_x, mu_y, g](real ay, real v) -> real {
    // Simple friction circle upper bound
    real ay_norm = ay / g;
    return g * std::sqrt(-ay_norm*ay_norm + mu_y*mu_y)*mu_x/mu_y;
  };

  auto GG_shape2 = [mu_x, mu_y, g](real ay, real v) -> real {
    // Simple friction circle lower bound
    real ay_norm = ay / g;
    return -g * std::sqrt(-ay_norm*ay_norm + mu_y*mu_y)*mu_x/mu_y;
  };

  auto gg_range_min = [mu_y, g](real v) -> real {
    return -mu_y * g;
  };

  auto gg_range_max = [mu_y, g](real v) -> real {
    return +mu_y * g;
  };

  GgRangeMaxMin gg_range = {gg_range_min, gg_range_max};

  // 4. Create Fb2d object
  Fb2d fbga(GG_shape1, GG_shape2, gg_range);

  // 5. Compute optimal velocity profile
  real total_time = fbga.compute(SS_eval_vec, KK_eval_vec, v_initial);

  // 6. Extract results
  std::vector<real> VX_eval_vec(SS_eval_vec.size());
  std::vector<real> AX_eval_vec(SS_eval_vec.size());
  std::vector<real> AY_eval_vec(SS_eval_vec.size());

  for (size_t i = 0; i < SS_eval_vec.size(); ++i) {
    VX_eval_vec[i] = fbga.evalV(SS_eval_vec[i]);   // Velocity [m/s]
    AX_eval_vec[i] = fbga.evalAx(SS_eval_vec[i]);  // Longitudinal acceleration [m/s²]
    AY_eval_vec[i] = fbga.evalAy(SS_eval_vec[i]);  // Lateral acceleration [m/s²]
  }

  std::cout << "Total lap time: " << total_time << " seconds" << std::endl;

  std::cout << "Results:\n";
  std::cout << "SS: ";
  for (const auto& s : SS_eval_vec) {
    std::cout << s << " ";
  }
  std::cout << "\nKK: ";
  for (const auto& k : KK_eval_vec) {
    std::cout << k << " ";
  }
  std::cout << "\nVX: ";
  for (const auto& v : VX_eval_vec) {
    std::cout << v << " ";
  }
  std::cout << "\nAX: ";
  for (const auto& a : AX_eval_vec) {
    std::cout << a << " ";
  }
  std::cout << "\nAY: ";
  for (const auto& a : AY_eval_vec) {
    std::cout << a << " ";
  }
  std::cout << std::endl;
  return 0;
}
