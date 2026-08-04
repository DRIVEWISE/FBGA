# FBGA Usage Guide

This guide covers both forward-backward solvers in this repository:

- [Part 1: Fb2d (2D)](#part-1-fb2d-2d) -- `fb::fbga::Fb2d`, the 2D solver (curvature-only path, single G-G envelope).
- [Part 2: Fbga3dSolver (3D)](#part-2-fbga3dsolver-3d) -- `fb::fbga3d::Fbga3dSolver`, the 3D solver (banked/pitched paths, load transfer, adherence scaling), generic over the vehicle's GG-diagram model via `FbgaIndy`/`FbgaMoto`. Only built when configuring with `-DFBGA_BUILD_3D=ON` (see [README.md](README.md)).

---

## Part 1: Fb2d (2D)

### Overview

The `fb::fbga::Fb2d` (Forward-Backward) class is a sophisticated optimization tool for computing optimal velocity profiles along a race track or path, subject to physical constraints like tire grip, aerodynamic forces, and vehicle dynamics. It implements a forward-backward algorithm that finds the maximum velocity profile while respecting acceleration limits.

`Fb2d` lives in the `fbga` module (`fb::fbga` namespace); the `real`/`integer` vocabulary types and helper classes it builds on live in the `utils` module (`fb::utils` namespace). All the snippets below assume:

```cpp
using namespace fb::utils;
using namespace fb::fbga;
```

### Table of Contents

1. [Basic Concepts](#basic-concepts)
2. [Required Data Structures](#required-data-structures)
3. [Constructor Parameters](#constructor-parameters)
4. [Basic Usage Example](#basic-usage-example)
5. [Advanced Configuration](#advanced-configuration)
6. [API Reference](#api-reference)
7. [Complete Example](#complete-example)

### Basic Concepts

The algorithm requires:

- A path discretized into segments with curvature information (curvilinear abscissa and curvature values).
- Upper and lower acceleration limit functions (G-G diagram)
- Lateral acceleration range constraints
- Initial velocity

### Required Data Structures

#### 1. Track Data or Race line

```cpp
std::vector<real> SS_vec;  // Arc length coordinates [m]
std::vector<real> KK_vec;  // Curvature values [1/m]
real v_initial;            // Initial velocity [m/s]
```

#### 2. Constraint Functions

You need to define three key functions:

##### Upper Acceleration Bound

```cpp
auto GG_shape1 = [](real ay, real v) -> real {
    // Return maximum longitudinal acceleration [m/s²]
    // given lateral acceleration ay [m/s²] and velocity v [m/s]
    return max_ax_function(ay, v);
};
```

##### Lower Acceleration Bound

```cpp
auto GG_shape2 = [](real ay, real v) -> real {
    // Return minimum (most negative) longitudinal acceleration [m/s²]
    // given lateral acceleration ay [m/s²] and velocity v [m/s]
    return min_ax_function(ay, v);
};
```

##### Lateral Acceleration Range

```cpp
auto gg_range_min = [](real v) -> real {
    // Return minimum lateral acceleration [m/s²] at velocity v [m/s]
    return min_ay_function(v);
};

auto gg_range_max = [](real v) -> real {
    // Return maximum lateral acceleration [m/s²] at velocity v [m/s]
    return max_ay_function(v);
};

GgRangeMaxMin gg_range = {gg_range_min, gg_range_max};
```

### Constructor Parameters

```cpp
Fb2d fb2d(
    const std::function<real(real, real)> &gg_Upper,  // Upper acceleration bound
    const std::function<real(real, real)> &gg_Lower,  // Lower acceleration bound
    const GgRangeMaxMin &gg_range                      // Lateral acceleration range
);
```

### Basic Usage Example

```cpp
#include <fbga/fb2d.hh>
#include <utils/types.hh>

using namespace fb::utils;
using namespace fb::fbga;

int main() {
    // 1. Define track data
    std::vector<real> SS_vec = {0.0, 5.0, 10.0, 15.0, 20.0};  // Arc length [m]
    std::vector<real> KK_vec = {0.0, 0.001, 0.05, 0.04, 0.0};    // Curvature [1/m]
    real v_initial = 20.0;  // Initial velocity [m/s]
    
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
    Fb2d fb2d(GG_shape1, GG_shape2, gg_range);
    
    // 5. Compute optimal velocity profile
    real total_time = fb2d.compute(SS_vec, KK_vec, v_initial);
    
    // 6. Extract results
    std::vector<real> VX_vec(SS_vec.size());
    std::vector<real> AX_vec(SS_vec.size());
    std::vector<real> AY_vec(SS_vec.size());
    
    for (size_t i = 0; i < SS_vec.size(); ++i) {
        VX_vec[i] = fb2d.evalV(SS_vec[i]);   // Velocity [m/s]
        AX_vec[i] = fb2d.evalAx(SS_vec[i]);  // Longitudinal acceleration [m/s²]
        AY_vec[i] = fb2d.evalAy(SS_vec[i]);  // Lateral acceleration [m/s²]
    }
    
    std::cout << "Total lap time: " << total_time << " seconds" << std::endl;
    
    return 0;
}
```

### Advanced Configuration

#### Custom Solver Parameters

```cpp
// Define custom solver parameters (optional)
SolverParams custom_params;
custom_params.tolerance = 1.0e-12;      // Higher precision
custom_params.max_iter = 500;           // More iterations
custom_params.verbosity = "detailed";   // More verbose output

// Note: Currently, you need to modify the Fb2d class to accept custom solver parameters
// or set them after construction if the API allows it
```

#### Vehicle-Specific Constraint Functions

For a motorcycle (as in the example), you might have more complex constraints:

```cpp
struct VehicleParams {
    real b = 0.73;       // Wheelbase [m]
    real h = 0.69;       // Center of mass height [m]
    real mu_X = 1.30;    // Longitudinal friction
    real mu_Y = 1.40;    // Lateral friction
    real m = 250.0;      // Mass [kg]
    real P = 145000.0;   // Maximum power [W]
    // ... aerodynamic coefficients, etc.
};

VehicleParams vehicle;

auto engine_limit = [&vehicle](real v) -> real {
    return vehicle.P / (vehicle.m * v);  // Power-limited acceleration
};

auto GG_shape1 = [&vehicle, &engine_limit](real ay, real v) -> real {
    real ax_engine = engine_limit(v);
    real ax_friction = friction_limit(ay, v, vehicle);
    real ax_wheeling = wheeling_limit(ay, v, vehicle);
    
    return std::min({ax_engine, ax_friction, ax_wheeling});
};
```

### API Reference

#### Core Methods

##### `compute()`
```cpp
real compute(std::vector<real> const &SS, std::vector<real> const &KK, real v0, real vfmax = 130.0);
```
- **Parameters:**
  - `SS`: Arc length vector [m]
  - `KK`: Curvature vector [1/m]
  - `v0`: Initial velocity [m/s]
  - `vfmax`: Maximum allowed final velocity [m/s] (optional, defaults to the class's `VMAX_SPEED` constant, currently 130.0)
- **Returns:** Total time [s]
- **Description:** Main computation method that runs the forward-backward algorithm

##### Evaluation Methods
```cpp
real evalV(real s) const;      // Velocity at arc length s [m/s]
real evalAx(real s) const;     // Longitudinal acceleration at s [m/s²]
real evalAy(real s) const;     // Lateral acceleration at s [m/s²]
real evalVmax(real s) const;   // Maximum velocity at s [m/s]
```

### Tips and Best Practices

1. **Mesh Resolution**: Use appropriate spacing in your `SS_vec`. Too coarse and you'll miss important features; too fine and computation becomes slow.

2. **Initial Velocity**: Choose a reasonable initial velocity. Too high and the algorithm may fail to find a solution.

3. **Constraint Functions**: Make sure your upper and lower bounds are consistent (upper ≥ lower for all valid inputs).

4. **Debugging**: Use `get_dump()` method to identify problematic segments if the algorithm fails to converge.

5. **Performance**: The algorithm complexity scales with the number of segments, so balance accuracy with computational cost.

### Common Issues and Solutions

#### Algorithm Doesn't Converge

- Check that constraint functions are well-defined for all input ranges
- Verify that the track data (curvature) is reasonable
- Reduce tolerance or increase maximum iterations of the solver
- Check for discontinuities in constraint functions

#### Unrealistic Results

- Verify units (make sure everything is in SI: meters, seconds, m/s², etc.)
- Check constraint function implementations
- Validate input track data for physical reasonableness

#### Performance Issues

- Reduce mesh density if acceptable for your application
- Optimize constraint function calculations
- Consider using simpler constraint models for initial testing

---

## Part 2: Fbga3dSolver (3D)

### Overview

`fb::fbga3d::Fbga3dSolver<GggvModel>` is the 3D counterpart to `Fb2d`: a forward-backward velocity-profile solver for paths with banking, pitch, and load transfer, not just curvature. It's a class template generic over the vehicle's GG-diagram model (`GggvModel`); this repository ships two instantiations as type aliases:

```cpp
using FbgaIndy = Fbga3dSolver<GggvIndy>;  // spline-based envelope, loaded from ./data/INDY/*.npy
using FbgaMoto = Fbga3dSolver<GggvMoto>;  // closed-form aero/adherence envelope, no external data
```

Both share the exact same public API -- only the GG-diagram model differs. `fbga3d` lives in the `fb::fbga3d` namespace and depends on `fb::utils`/`fb::solvers` like `fbga` does. All the snippets below assume:

```cpp
using namespace fb::utils;
using namespace fb::fbga3d;
```

Only built when configuring CMake with `-DFBGA_BUILD_3D=ON` (see [README.md](README.md)); ported from and numerically verified against the original [FBGA_3D](https://github.com/DRIVEWISE/FBGA_3D) repository (see `.claude/FBGA3D_INTEGRATION_PLAN.md` for the porting/audit history). Plotting (`eval_shell_plot*`) and constraint-satisfaction-checking helpers from the original FBGA_3D are not ported -- this is the core solve path only.

### Basic Concepts

The algorithm requires:

- A road/path discretized into nodes with, at each node: arc-length abscissa, lateral offset `n` and heading offset `chi` from a reference line, and the reference line's own pitch/roll/yaw angles (`mu`/`phi`/`theta`) and their geometric derivatives.
- A vehicle's GG-diagram model (`GggvIndy` or `GggvMoto`, or your own type implementing the same interface -- see below).
- Per-node adherence scaling (`alpha`, typically `1.0`).
- Initial velocity.

### Required Data Structures

#### 1. Trajectory Data

```cpp
TrajectoryOffsetAndAnglesContainer TOA;
TOA.reference.abscissa    = {...};  // Arc length [m]
TOA.reference.mu          = {...};  // Pitch angle [rad]
TOA.reference.phi         = {...};  // Roll angle [rad]
TOA.reference.theta       = {...};  // Yaw angle [rad]
TOA.reference.mu_prime    = {...};  // d(mu)/ds
TOA.reference.phi_prime   = {...};  // d(phi)/ds
TOA.reference.theta_prime = {...};  // d(theta)/ds (curvature-equivalent term)
TOA.offset.n              = {...};  // Lateral offset from the reference line [m]
TOA.offset.chi            = {...};  // Heading offset from the reference line [rad]
TOA.adherence.alpha       = {...};  // Adherence scaling, typically 1.0
```

All vectors must be the same length (one entry per node; a track with `N` nodes has `N - 1` solved cells).

#### 2. GG-Diagram Model (`GggvModel`)

`GggvIndy` loads a bilinear-interpolated envelope from `.npy` files (relative to the working directory: `./data/INDY/g_list.npy`, `v_list.npy`, `ay_max.npy`, `ax_min.npy`, `ax_max.npy`) plus an engine-power curve:

```cpp
GggvIndy gggv;                                    // default: reads ./data/INDY/*.npy
GggvIndy gggv(scaling_factors);                    // + custom ScalingGggvFactors
GggvIndy gggv(spline_data_collection);              // fully custom spline data, no file I/O
```

`GggvMoto` is closed-form (wheeling/stoppie/adherence/aero constraints from `MotoData`), no external data needed:

```cpp
GggvMoto gggv;  // default MotoData
```

Any type implementing `a_x_push`/`a_x_pull`/`a_x_neutral`/`a_x_eng`/`a_y_lim` with the same signatures as `GggvIndy`/`GggvMoto` can be used as the `GggvModel` template parameter for a new `Fbga3dSolver<YourModel>` instantiation.

### Constructor and Setup

```cpp
FbgaIndy solver;                    // default-constructs its own GggvIndy
solver.model() = GggvIndy(spline_data);  // or: replace/configure the model directly
solver.set_max_velocity(90.0);      // optional cap, default 130.0 m/s
```

### Basic Usage Example

```cpp
#include <fbga3d/fbga_moto.hh>

using namespace fb::utils;
using namespace fb::fbga3d;

int main() {
    // 1. Define trajectory data (flat, non-banked track with a simple curvature profile)
    TrajectoryOffsetAndAnglesContainer TOA;
    TOA.reference.abscissa    = {0.0, 50.0, 100.0, 150.0, 200.0};
    TOA.reference.mu          = {0.0, 0.0, 0.0, 0.0, 0.0};
    TOA.reference.phi         = {0.0, 0.0, 0.0, 0.0, 0.0};
    TOA.reference.theta       = {0.0, 0.05, 0.55, 0.8, 0.8};
    TOA.reference.mu_prime    = {0.0, 0.0, 0.0, 0.0, 0.0};
    TOA.reference.phi_prime   = {0.0, 0.0, 0.0, 0.0, 0.0};
    TOA.reference.theta_prime = {0.001, 0.01, 0.005, -0.01, 0.0};
    TOA.offset.n   = {0.0, 0.0, 0.0, 0.0, 0.0};
    TOA.offset.chi = {0.0, 0.0, 0.0, 0.0, 0.0};
    TOA.adherence.alpha = {1.0, 1.0, 1.0, 1.0, 1.0};

    real v_initial = 20.0;  // Initial velocity [m/s]

    // 2. Create the solver (FbgaMoto needs no external data files)
    FbgaMoto solver;

    // 3. Compute optimal velocity profile
    real total_time = solver.compute(TOA, v_initial);

    // 4. Extract results
    for (real s : TOA.reference.abscissa) {
        real v = solver.eval_V(s);              // Velocity [m/s]
        real ax = solver.eval_A_tilde_x(s);      // Longitudinal accel (path frame) [m/s²]
        real ay = solver.eval_A_tilde_y(s);      // Lateral accel (path frame) [m/s²]
        (void) v; (void) ax; (void) ay;
    }

    std::cout << "Total lap time: " << total_time << " seconds" << std::endl;

    return 0;
}
```

See `examples/fbga3d/yas_marina_indy.cc` for a complete example driving `FbgaIndy` off real circuit CSV/`.npy` data, and `test/fbga3d/test_fbga_solver.cc` for the `FbgaMoto`/`FbgaIndy` scenario this example mirrors.

### Advanced Configuration

#### Yellow-Flag (Speed-Limited) Segments

```cpp
YellowFlagData yellow_flag;
yellow_flag.is_yellow  = true;
yellow_flag.v_des_max  = 30.0;   // Target speed under yellow flag [m/s]
yellow_flag.a_des_min  = -5.0;   // Target braking deceleration [m/s²]

real total_time = solver.compute(TOA, v_initial, yellow_flag);
```

### API Reference

#### Core Methods

##### `compute()`
```cpp
real compute(TrajectoryOffsetAndAnglesContainer const &TOA, real V0 = -1.0, YellowFlagData const &yellow_flag = YellowFlagData());
```
- **Parameters:**
  - `TOA`: Trajectory data (see above)
  - `V0`: Initial velocity [m/s] (optional; `<= 0` falls back to `max_velocity`)
  - `yellow_flag`: Yellow-flag data (optional, defaults to inactive)
- **Returns:** Total time [s]
- **Description:** Main computation method that runs the forward-backward algorithm (`BY()` if yellow-flag active, then `FW()`, then `BW()`)

##### `get_solution()`
```cpp
SolutionContainer get_solution() const;  // per-cell V0, V1, V_dot
```

##### `get_dump()`
```cpp
std::vector<int> get_dump() const;  // indices of cells that failed to converge (debugging aid)
```

##### Evaluation Methods (by arc-length `s`)
```cpp
real eval_V(real s) const;           // Velocity [m/s]
real eval_V_dot(real s) const;       // Velocity time-derivative [m/s²]
real eval_Vmax(real s) const;        // Maximum reachable velocity [m/s]
real eval_A_hat_x(real s) const;     // Longitudinal accel (vehicle frame) [m/s²]
real eval_A_hat_y(real s) const;     // Lateral accel (vehicle frame) [m/s²]
real eval_A_hat_z(real s) const;     // Vertical accel (vehicle frame) [m/s²]
real eval_A_tilde_x(real s) const;   // Longitudinal accel minus gravity component [m/s²]
real eval_A_tilde_y(real s) const;   // Lateral accel minus gravity component [m/s²]
real eval_A_tilde_z(real s) const;   // Vertical accel minus gravity component [m/s²]
real eval_alpha(real s) const;       // Interpolated adherence scaling
SegmentType eval_segment_type(real s) const;  // FORWARD/BACKWARD/TRANSITION/YELLOWFLAG/..._NAN
```

### Tips and Best Practices

1. **Mesh Resolution**: Same trade-off as `Fb2d` -- too coarse misses features, too fine is slow. `FbgaIndy`/`FbgaMoto` scale to several thousand nodes fine (a 5214-node Yas Marina lap solves in ~15-20ms).

2. **`GggvIndy` needs its data files**: run from a working directory where `./data/INDY/*.npy` resolves, or construct it from an explicit `SplineDataCollection` instead of the default constructor.

3. **Flat/non-banked tracks are a valid special case**: set `mu`/`phi`/`mu_prime`/`phi_prime` to zero and this reduces to a 2D-equivalent problem, useful for sanity-checking against `Fb2d` on the same curvature profile.

4. **Debugging**: use `get_dump()` the same way as `Fb2d` -- a cell appearing there had its forward/backward/transition solve fail to converge and fell back to a clamped value; check `eval_segment_type()` at that `s` for which pass (`FORWARD_NAN`/`BACKWARD_NAN`/`YELLOWFLAG_NAN`) failed.

### Common Issues and Solutions

#### `GggvIndy` throws or default-constructs into garbage data

- Check the working directory: `setup_std()` reads paths relative to the CWD (`./data/INDY/...`), not relative to the source file or executable.
- Prefer constructing from an explicit `SplineDataCollection` when the CWD isn't controlled by you (e.g. library embedded in another application).

#### Algorithm Doesn't Converge / Cells in `get_dump()`

- Same checklist as `Fb2d`: verify the GG-diagram model is well-defined (no NaNs) across the velocity/lateral-accel/vertical-accel ranges actually reached by the trajectory.
- Check that `mu`/`phi`/`theta` and their derivatives are consistent with `n`/`chi` (e.g. no extreme banking combined with a large lateral offset that pushes the vehicle off the physically valid range).
