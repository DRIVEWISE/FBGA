#include <iostream>

#include <fbga3d/gggv_indy.hh>
#include <fbga3d/gggv_moto.hh>
#include <utils/types.hh>

using namespace fb::fbga3d;

int main()
{
  std::cout << "FBGA3D - vehicle GG-diagram models\n";
  std::cout << " > Running example: gggv_demo\n";
  std::cout << " > Note: this demoes the ported GggvIndy/GggvMoto GG-diagram models only.\n";
  std::cout << " >       The FbgaIndy/FbgaMoto forward-backward solvers built on top of them\n";
  std::cout << " >       are not ported yet (see FBGA3D_INTEGRATION_PLAN.md).\n\n";

  // --- GggvMoto: closed-form aero/adherence model, no external data needed ---
  std::cout << "=== GggvMoto ===\n";
  GggvMoto gggv_moto;

  std::cout << "Engine-limited longitudinal acceleration vs speed:\n";
  for (real V = 10.0; V <= 50.0; V += 10.0)
  {
    std::cout << "  V: " << V << " m/s, a_x_eng: " << gggv_moto.a_x_eng(V) << " m/s^2\n";
  }

  {
    const real V = 20.0;
    const real az_tilde = 9.81;
    const real a_y_lim = gggv_moto.a_y_lim(V, az_tilde);
    std::cout << "\nGG envelope at V = " << V << " m/s, az_tilde = " << az_tilde << " m/s^2:\n";
    for (real ay = -a_y_lim; ay <= a_y_lim; ay += a_y_lim / 2.0)
    {
      std::cout << "  ay: " << ay
                << ", a_x_max: " << gggv_moto.a_x_push(ay, V, az_tilde)
                << ", a_x_min: " << gggv_moto.a_x_pull(ay, V, az_tilde) << "\n";
    }
  }

  // --- GggvIndy: bilinear-interpolated spline model, loads ./data/INDY/*.npy ---
  std::cout << "\n=== GggvIndy ===\n";
  std::cout << "(loading spline data from ./data/INDY/ -- run this example from the repo root)\n";
  GggvIndy gggv_indy;

  {
    const real V = 30.0;
    const real az_tilde = 10.0;
    const real a_y_lim = gggv_indy.a_y_lim(V, az_tilde);
    std::cout << "GG envelope at V = " << V << " m/s, az_tilde = " << az_tilde << " m/s^2:\n";
    for (real ay = -a_y_lim; ay <= a_y_lim; ay += a_y_lim / 2.0)
    {
      std::cout << "  ay: " << ay
                << ", a_x_max: " << gggv_indy.a_x_push(ay, V, az_tilde)
                << ", a_x_min: " << gggv_indy.a_x_pull(ay, V, az_tilde) << "\n";
    }
  }

  return 0;
}
