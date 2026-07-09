#!/usr/bin/env python3
"""
Example: use the FWBW solver from Python with lambda-defined G-G bounds,
and plot the resulting velocity/acceleration profiles and the G-G diagram.

Scenario mirrors src_tests/GIGI_test_basic.cc (friction-circle bounds).
"""

import math

import matplotlib.pyplot as plt
import numpy as np

import fbga_py as fb


def main():
    print("FBGA - Forward Backward Generic Acceleration constraints")
    print(" > Running FWBW example with Python lambda functions")

    # 1. Define track data
    SS_vec = [0.0, 50.0, 100.0, 150.0, 200.0]           # Arc length [m]
    KK_vec = [0.001, 0.01, 0.005, -0.01, 0.0]           # Curvature [1/m]
    v_initial = 20.0                                     # Initial velocity [m/s]

    print(" > Track data:")
    print(f" >   SS_vec: {len(SS_vec)} points")
    print(f" >   KK_vec: {len(KK_vec)} points")
    print(f" >   v_initial: {v_initial} m/s")

    # Interpolate track data to get evaluation points
    numpt_eval = math.ceil(SS_vec[-1])
    SS_eval_vec = [float(i) for i in range(numpt_eval)]

    KK_eval_vec = []
    for s in SS_eval_vec:
        for j in range(len(SS_vec) - 1):
            if SS_vec[j] <= s <= SS_vec[j + 1]:
                t = (s - SS_vec[j]) / (SS_vec[j + 1] - SS_vec[j])
                KK_eval_vec.append(KK_vec[j] + t * (KK_vec[j + 1] - KK_vec[j]))
                break
        else:
            KK_eval_vec.append(KK_vec[-1])

    # 2. Define vehicle parameters (friction circle)
    mu_x = 1.3   # Longitudinal friction coefficient
    mu_y = 1.4   # Lateral friction coefficient
    g = fb.GRAVITY

    print(" > Defining G-G diagram functions with Python lambdas...")

    def gg_upper(ay: float, v: float) -> float:
        """Friction-circle upper bound on ax."""
        ay_norm = ay / g
        return g * math.sqrt(max(0.0, mu_y * mu_y - ay_norm * ay_norm)) * mu_x / mu_y

    def gg_lower(ay: float, v: float) -> float:
        """Friction-circle lower bound on ax."""
        return -gg_upper(ay, v)

    def ay_max_func(v: float) -> float:
        return mu_y * g

    def ay_min_func(v: float) -> float:
        return -mu_y * g

    # Create the range struct
    gg_range = fb.GGRangeMaxMin()
    gg_range.min = ay_min_func
    gg_range.max = ay_max_func

    # 3. Create FWBW object with lambda functions
    print(" > Creating FWBW instance with lambda functions...")
    fwbw_solver = fb.FWBW(gg_upper, gg_lower, gg_range)

    # 4. Test the functions by evaluating at test points
    test_v, test_ay = 10.0, 2.0
    print(f"   At v={test_v} m/s, ay={test_ay} m/s^2:")
    print(f"   Max ax: {gg_upper(test_ay, test_v):.3f} m/s^2")
    print(f"   Min ax: {gg_lower(test_ay, test_v):.3f} m/s^2")
    print(f"   Max ay: {ay_max_func(test_v):.3f} m/s^2")
    print(f"   Min ay: {ay_min_func(test_v):.3f} m/s^2")

    # 5. Compute optimal velocity profile
    print(f" > Computing optimal trajectory for {len(SS_eval_vec)} points...")
    print(f" > Track length: {SS_vec[-1]} m")
    print(f" > Initial velocity: {v_initial} m/s")

    total_time = fwbw_solver.compute(SS_eval_vec, KK_eval_vec, v_initial)

    # 6. Extract results
    AX_eval_vec, AY_eval_vec, VX_eval_vec = fwbw_solver.evaluate(SS_eval_vec)

    print(f" > Total lap time: {total_time:.3f} seconds")
    print(f" > Final velocity: {VX_eval_vec[-1]:.3f} m/s")
    print(f" > Max velocity reached: {max(VX_eval_vec):.3f} m/s")
    print(" > Optimization completed successfully!")

    # 7. Create visualizations
    create_plots(
        SS_eval_vec, VX_eval_vec, AX_eval_vec, AY_eval_vec, KK_eval_vec,
        gg_upper, gg_lower, ay_max_func, ay_min_func,
    )


def create_plots(SS_eval_vec, VX_eval_vec, AX_eval_vec, AY_eval_vec, KK_eval_vec,
                  gg_upper, gg_lower, ay_max_func, ay_min_func):
    """Create comprehensive plots for the FWBW results."""

    print("\nCreating plots...")

    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 10))

    # Plot 1: Velocity vs Position
    ax1.plot(SS_eval_vec, VX_eval_vec, "b-", linewidth=2, label="Velocity")
    ax1.set_xlabel("Position s [m]")
    ax1.set_ylabel("Velocity [m/s]")
    ax1.set_title("Optimal Velocity Profile")
    ax1.grid(True, alpha=0.3)
    ax1.legend()

    # Plot 2: Longitudinal Acceleration vs Position
    ax2.plot(SS_eval_vec, AX_eval_vec, "r-", linewidth=2, label="Longitudinal Acceleration")
    ax2.set_xlabel("Position s [m]")
    ax2.set_ylabel("Longitudinal Acceleration [m/s^2]")
    ax2.set_title("Longitudinal Acceleration Profile")
    ax2.grid(True, alpha=0.3)
    ax2.legend()
    ax2.axhline(y=0, color="k", linestyle="--", alpha=0.5)

    # Plot 3: Lateral Acceleration vs Position (with curvature on twin axis)
    ax3.plot(SS_eval_vec, AY_eval_vec, "g-", linewidth=2, label="Lateral Acceleration")
    ax3.set_xlabel("Position s [m]")
    ax3.set_ylabel("Lateral Acceleration [m/s^2]")
    ax3.set_title("Lateral Acceleration Profile")
    ax3.grid(True, alpha=0.3)
    ax3.legend(loc="upper left")
    ax3.axhline(y=0, color="k", linestyle="--", alpha=0.5)

    ax3_twin = ax3.twinx()
    ax3_twin.plot(SS_eval_vec, KK_eval_vec, "orange", linestyle=":", alpha=0.7, label="Curvature")
    ax3_twin.set_ylabel("Curvature [1/m]", color="orange")
    ax3_twin.tick_params(axis="y", labelcolor="orange")
    ax3_twin.legend(loc="upper right")

    # Plot 4: G-G Diagram
    create_gg_diagram(ax4, AX_eval_vec, AY_eval_vec, VX_eval_vec,
                       gg_upper, gg_lower, ay_max_func, ay_min_func)

    plt.tight_layout()
    out_path = "fwbw_example_plots.png"
    plt.savefig(out_path, dpi=150)
    print(f" > Saved 2D plots to {out_path}")
    plt.show()

    # 3D plot in (ay, v, ax) space
    create_3d_plot(AY_eval_vec, VX_eval_vec, AX_eval_vec,
                   gg_upper, gg_lower, ay_max_func, ay_min_func)


def create_gg_diagram(ax, AX_eval_vec, AY_eval_vec, VX_eval_vec,
                       gg_upper, gg_lower, ay_max_func, ay_min_func):
    """Create G-G diagram showing trajectory and boundaries at a few velocities."""

    scatter = ax.scatter(AX_eval_vec, AY_eval_vec, c=VX_eval_vec, cmap="viridis",
                          s=20, alpha=0.7, label="Trajectory")
    ax.plot(AX_eval_vec[0], AY_eval_vec[0], "go", markersize=8, label="Start")
    ax.plot(AX_eval_vec[-1], AY_eval_vec[-1], "ro", markersize=8, label="End")

    velocities = [5, 10, 15, 20]
    colors = ["lightblue", "lightgreen", "lightyellow", "lightcoral"]

    for v, color in zip(velocities, colors):
        ay_range = np.linspace(ay_min_func(v), ay_max_func(v), 100)
        ax_upper = [gg_upper(ay, v) for ay in ay_range]
        ax_lower = [gg_lower(ay, v) for ay in ay_range]

        ay_envelope = np.concatenate([ay_range, ay_range[::-1]])
        ax_envelope = np.concatenate([ax_upper, ax_lower[::-1]])

        ax.fill(ax_envelope, ay_envelope, alpha=0.2, color=color,
                label=f"G-G envelope v={v} m/s")

    ax.set_xlabel("Longitudinal Acceleration [m/s^2]")
    ax.set_ylabel("Lateral Acceleration [m/s^2]")
    ax.set_title("G-G Diagram with Lambda-defined Boundaries")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=7)
    ax.axis("equal")

    plt.colorbar(scatter, ax=ax, label="Velocity [m/s]")

    ax.axhline(y=0, color="k", linestyle="-", alpha=0.3)
    ax.axvline(x=0, color="k", linestyle="-", alpha=0.3)


def create_3d_plot(AY_eval_vec, VX_eval_vec, AX_eval_vec,
                    gg_upper, gg_lower, ay_max_func, ay_min_func):
    """Create 3D plot showing the trajectory in (ay, v, ax) space with G-G bound surfaces."""

    fig = plt.figure(figsize=(12, 8))
    ax = fig.add_subplot(111, projection="3d")

    ax.plot(AY_eval_vec, VX_eval_vec, AX_eval_vec, "bo-", markersize=3,
            linewidth=2, label="Trajectory", alpha=0.8)
    ax.scatter(AY_eval_vec[0], VX_eval_vec[0], AX_eval_vec[0],
               color="green", s=100, label="Start")
    ax.scatter(AY_eval_vec[-1], VX_eval_vec[-1], AX_eval_vec[-1],
               color="red", s=100, label="End")

    v_range = np.linspace(2, 20, 20)
    ay_samples = 50

    for v in v_range[::3]:
        ay_range = np.linspace(ay_min_func(v), ay_max_func(v), ay_samples)
        ax_upper = [gg_upper(ay, v) for ay in ay_range]
        ax_lower = [gg_lower(ay, v) for ay in ay_range]

        ax.plot(ay_range, [v] * len(ay_range), ax_upper, "r-", alpha=0.3)
        ax.plot(ay_range, [v] * len(ay_range), ax_lower, "b-", alpha=0.3)

    ax.set_xlabel("Lateral Acceleration ay [m/s^2]")
    ax.set_ylabel("Velocity v [m/s]")
    ax.set_zlabel("Longitudinal Acceleration ax [m/s^2]")
    ax.set_title("3D G-G-V Space with Lambda-defined Boundaries")
    ax.legend()
    ax.view_init(elev=20, azim=45)

    plt.tight_layout()
    out_path = "fwbw_example_3d.png"
    plt.savefig(out_path, dpi=150)
    print(f" > Saved 3D plot to {out_path}")
    plt.show()


if __name__ == "__main__":
    try:
        main()
        print("\nFWBW example completed successfully!")
    except Exception as exc:
        print(f"Error: {exc}")
        print("Make sure the fbga_py module is built and installed correctly.")
        import traceback
        traceback.print_exc()
