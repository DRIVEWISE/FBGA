#include <fbga3d/fbga_indy.hh>

#include <utils/gg_utils.hh>

#include <cxxopts.hpp>
#include <rapidcsv.h>
#include <npy.hpp>

#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

// Reproduces FBGA_3D's FBGA3D_test_INDY.cc scenario against the ported FbgaIndy, dumping the
// same per-node columns so they can be diffed against a reference run of the original FBGA_INDY
// on the same data/options. Defaults mirror FBGA3D_test_INDY.cc's own.
using namespace fb::fbga3d;

int main(int argc, char *argv[])
{
  cxxopts::Options options("fbga3d_audit_indy_yas_marina");
  options.add_options()
    ("c,circuit", "Circuit CSV file name", cxxopts::value<std::string>()->default_value("Yas_Marina_raceline.csv"))
    ("f,folder", "Folder path of the circuit data", cxxopts::value<std::string>()->default_value("./data/INDY/"))
    ("y,yellow", "Yellow flag active", cxxopts::value<bool>()->default_value("false"))
    ("s,speed", "Target speed, yellow flag [m/s]", cxxopts::value<real>()->default_value("100"))
    ("a,acceleration", "Target acceleration, yellow flag [m/s^2]", cxxopts::value<real>()->default_value("-5"))
    ("v,initial-speed", "Initial speed [m/s]", cxxopts::value<real>()->default_value("20"))
    ("b,bound", "Use the wider engine envelope / scaling factors", cxxopts::value<bool>()->default_value("false"))
    ("o,output", "Write per-node output CSV", cxxopts::value<bool>()->default_value("false"))
    ("n,note", "Suffix appended to the output file name", cxxopts::value<std::string>()->default_value("ported"))
    ("h,help", "Print help");
  const auto result = options.parse(argc, argv);
  if (result.count("help"))
  {
    std::cout << options.help() << "\n";
    return 0;
  }

  const std::string circuit_name = result["circuit"].as<std::string>();
  const std::string folder_path = result["folder"].as<std::string>();
  const bool yellow_flag_active = result["yellow"].as<bool>();
  const real v_des_yf = result["speed"].as<real>();
  const real a_des_yf = result["acceleration"].as<real>();
  const real v_initial = result["initial-speed"].as<real>();
  const bool bound = result["bound"].as<bool>();
  const bool output_write = result["output"].as<bool>();
  const std::string note = result["note"].as<std::string>();

  const std::string road_file_path = folder_path + circuit_name;

  rapidcsv::Document doc(
    road_file_path,
    rapidcsv::LabelParams(0, -1),
    rapidcsv::SeparatorParams(',', true),
    rapidcsv::ConverterParams(),
    rapidcsv::LineReaderParams(true, '#')
  );

  const std::vector<std::string> column_names = {
    "n_rl_m", "chi_rl_rad", "s_ref_rl_m", "theta_ref_rl_rad", "mu_ref_rl_rad",
    "phi_ref_rl_rad", "dtheta_ref_rl_radpm", "dmu_ref_rl_radpm", "dphi_ref_rl_radpm"
  };

  std::vector<std::vector<real>> YVEC_SPLNE;
  YVEC_SPLNE.reserve(column_names.size());
  for (const auto &col_name : column_names)
  {
    YVEC_SPLNE.push_back(doc.GetColumn<real>(col_name));
  }

  const std::vector<real> Sref = doc.GetColumn<real>("s_ref_rl_m");

  const fb::utils::LinearInterpolatorSet traj_spline(Sref, YVEC_SPLNE, column_names);

  const std::vector<real> g_data = std::move(npy::read_npy<real>(folder_path + "g_list.npy").data);
  const std::vector<real> v_data = std::move(npy::read_npy<real>(folder_path + "v_list.npy").data);
  const std::vector<real> ay_max_data = std::move(npy::read_npy<real>(folder_path + "ay_max.npy").data);
  const std::vector<real> ax_min_data = std::move(npy::read_npy<real>(folder_path + "ax_min.npy").data);
  const std::vector<real> ax_max_data = std::move(npy::read_npy<real>(folder_path + "ax_max.npy").data);

  // Matches FBGA3D_test_INDY.cc's own -b/--bound switch between the two hand-tuned envelopes.
  const std::vector<real> v_eng_data = {0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0};
  std::vector<real> ax_eng_data;
  ScalingGggvFactors scaling_factors;
  if (bound)
  {
    ax_eng_data = {10.0, 10.0, 9.5, 9.0, 7.4, 4.9, 2.9, 0.6, 0.1, 0.0};
    scaling_factors = {1.2, 1.35, 1.1, 1.55, 1.80};
  }
  else
  {
    ax_eng_data = {5.0, 10.0, 13.0, 10.0, 7.0, 4.6, 2.5, 0.3, 0.0, 0.0};
    scaling_factors = {1.0, 1.0, 1.0, 1.3, 1.3};
  }

  const SplineDataCollection spline_data{
    v_data, g_data, ay_max_data, ax_min_data, ax_max_data, scaling_factors, v_eng_data, ax_eng_data
  };

  GggvIndy gggv_indy(spline_data);
  FbgaIndy fbga_indy;
  fbga_indy.model() = gggv_indy;

  const integer numpts_gg = static_cast<integer>(Sref.size());
  std::vector<real> S_gg(numpts_gg);
  std::vector<real> THETA_gg(numpts_gg), PHI_GG(numpts_gg), MU_GG(numpts_gg);
  std::vector<real> dTHETA_gg(numpts_gg), dPHI_GG(numpts_gg), dMU_GG(numpts_gg);
  std::vector<real> CHI_GG(numpts_gg), N_GG(numpts_gg);
  const std::vector<real> ALPHAS(static_cast<size_t>(numpts_gg), 1.0);

  for (integer i = 0; i < numpts_gg; ++i)
  {
    const real s = Sref[static_cast<size_t>(i)];
    S_gg[static_cast<size_t>(i)] = s;
    THETA_gg[static_cast<size_t>(i)] = traj_spline.eval("theta_ref_rl_rad", s);
    PHI_GG[static_cast<size_t>(i)] = traj_spline.eval("phi_ref_rl_rad", s);
    MU_GG[static_cast<size_t>(i)] = traj_spline.eval("mu_ref_rl_rad", s);
    dTHETA_gg[static_cast<size_t>(i)] = traj_spline.eval("dtheta_ref_rl_radpm", s);
    dPHI_GG[static_cast<size_t>(i)] = traj_spline.eval("dphi_ref_rl_radpm", s);
    dMU_GG[static_cast<size_t>(i)] = traj_spline.eval("dmu_ref_rl_radpm", s);
    CHI_GG[static_cast<size_t>(i)] = traj_spline.eval("chi_rl_rad", s);
    N_GG[static_cast<size_t>(i)] = traj_spline.eval("n_rl_m", s);
  }

  TrajectoryOffsetAndAnglesContainer TOA;
  TOA.offset.n = N_GG;
  TOA.offset.chi = CHI_GG;
  TOA.reference.mu = MU_GG;
  TOA.reference.phi = PHI_GG;
  TOA.reference.theta = THETA_gg;
  TOA.reference.mu_prime = dMU_GG;
  TOA.reference.phi_prime = dPHI_GG;
  TOA.reference.theta_prime = dTHETA_gg;
  TOA.reference.abscissa = S_gg;
  TOA.adherence.alpha = ALPHAS;

  YellowFlagData yellow_flag_data;
  yellow_flag_data.v_des_max = v_des_yf;
  yellow_flag_data.a_des_min = a_des_yf;
  yellow_flag_data.is_yellow = yellow_flag_active;

  const real T = fbga_indy.compute(TOA, v_initial, yellow_flag_data);

  std::cout << "Total lap time: " << std::setprecision(12) << T << " s\n";
  std::cout << "Num segments:   " << numpts_gg << "\n";

  if (!output_write)
  {
    return 0;
  }

  std::string base_circuit_name = circuit_name;
  const size_t last_dot = base_circuit_name.find_last_of('.');
  if (last_dot != std::string::npos)
  {
    base_circuit_name = base_circuit_name.substr(0, last_dot);
  }
  const std::string output_file_path = folder_path + "FBGA_" + base_circuit_name + "_output_" + note + ".csv";
  std::cout << "Writing output to " << output_file_path << "\n";
  std::ofstream output_file(output_file_path);
  output_file << "n_rl_m,chi_rl_rad,s_ref_rl_m,"
              << "theta_ref_rl_rad, mu_ref_rl_rad, phi_ref_rl_rad, dtheta_ref_rl_radpm, dmu_ref_rl_radpm, dphi_ref_rl_radpm,"
              << "alpha_lu,s_fb_m,v_fb_mps,axhat_fb_mps2,ayhat_fb_mps2,azhat_fb_mps2,axtilde_fb_mps2,aytilde_fb_mps2,aztilde_fb_mps2,v_dot_fb_mps2,v_max_fb_mps,g_x_fb_mps2,g_y_fb_mps2,g_z_fb_mps2,Omega_x_fb_radpm,Omega_y_fb_radpm,Omega_z_fb_radpm\n";
  output_file << std::fixed << std::setprecision(12);
  for (integer i = 0; i < numpts_gg; ++i)
  {
    const real s = S_gg[static_cast<size_t>(i)];
    output_file
      << traj_spline.eval("n_rl_m", s) << ","
      << traj_spline.eval("chi_rl_rad", s) << ","
      << traj_spline.eval("s_ref_rl_m", s) << ","
      << traj_spline.eval("theta_ref_rl_rad", s) << ","
      << traj_spline.eval("mu_ref_rl_rad", s) << ","
      << traj_spline.eval("phi_ref_rl_rad", s) << ","
      << traj_spline.eval("dtheta_ref_rl_radpm", s) << ","
      << traj_spline.eval("dmu_ref_rl_radpm", s) << ","
      << traj_spline.eval("dphi_ref_rl_radpm", s) << ","
      << fbga_indy.eval_alpha(s) << ","
      << s << ","
      << fbga_indy.eval_V(s) << ","
      << fbga_indy.eval_A_hat_x(s) << ","
      << fbga_indy.eval_A_hat_y(s) << ","
      << fbga_indy.eval_A_hat_z(s) << ","
      << fbga_indy.eval_A_tilde_x(s) << ","
      << fbga_indy.eval_A_tilde_y(s) << ","
      << fbga_indy.eval_A_tilde_z(s) << ","
      << fbga_indy.eval_V_dot(s) << ","
      << fbga_indy.eval_Vmax(s) << ","
      << fbga_indy.eval_g_x(s) << ","
      << fbga_indy.eval_g_y(s) << ","
      << fbga_indy.eval_g_z(s) << ","
      << fbga_indy.eval_Omega_x(s) << ","
      << fbga_indy.eval_Omega_y(s) << ","
      << fbga_indy.eval_Omega_z(s) << "\n";
  }

  return 0;
}
