# fb::fbga::Fb2d Usage Guide

## Overview

The `fb::fbga::Fb2d` (Forward-Backward) class is a sophisticated optimization tool for computing optimal velocity profiles along a race track or path, subject to physical constraints like tire grip, aerodynamic forces, and vehicle dynamics. It implements a forward-backward algorithm that finds the maximum velocity profile while respecting acceleration limits.

`Fb2d` lives in the `fbga` module (`fb::fbga` namespace); the `real`/`integer` vocabulary types and helper classes it builds on live in the `utils` module (`fb::utils` namespace). All the snippets below assume:

```cpp
using namespace fb::utils;
using namespace fb::fbga;
```

## Table of Contents

1. [Basic Concepts](#basic-concepts)
2. [Required Data Structures](#required-data-structures)
3. [Constructor Parameters](#constructor-parameters)
4. [Basic Usage Example](#basic-usage-example)
5. [Advanced Configuration](#advanced-configuration)
6. [API Reference](#api-reference)
7. [Complete Example](#complete-example)

## Basic Concepts

The algorithm requires:

- A path discretized into segments with curvature information (curvilinear abscissa and curvature values).
- Upper and lower acceleration limit functions (G-G diagram)
- Lateral acceleration range constraints
- Initial velocity

## Required Data Structures

### 1. Track Data or Race line

```cpp
std::vector<real> SS_vec;  // Arc length coordinates [m]
std::vector<real> KK_vec;  // Curvature values [1/m]
real v_initial;            // Initial velocity [m/s]
```

### 2. Constraint Functions

You need to define three key functions:

#### Upper Acceleration Bound

```cpp
auto GG_shape1 = [](real ay, real v) -> real {
    // Return maximum longitudinal acceleration [m/s²]
    // given lateral acceleration ay [m/s²] and velocity v [m/s]
    return max_ax_function(ay, v);
};
```

#### Lower Acceleration Bound

```cpp
auto GG_shape2 = [](real ay, real v) -> real {
    // Return minimum (most negative) longitudinal acceleration [m/s²]
    // given lateral acceleration ay [m/s²] and velocity v [m/s]
    return min_ax_function(ay, v);
};
```

#### Lateral Acceleration Range

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

## Constructor Parameters

```cpp
Fb2d fb2d(
    const std::function<real(real, real)> &gg_Upper,  // Upper acceleration bound
    const std::function<real(real, real)> &gg_Lower,  // Lower acceleration bound
    const GgRangeMaxMin &gg_range                      // Lateral acceleration range
);
```

## Basic Usage Example

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

## Advanced Configuration

### Custom Solver Parameters

```cpp
// Define custom solver parameters (optional)
SolverParams custom_params;
custom_params.tolerance = 1.0e-12;      // Higher precision
custom_params.max_iter = 500;           // More iterations
custom_params.verbosity = "detailed";   // More verbose output

// Note: Currently, you need to modify the Fb2d class to accept custom solver parameters
// or set them after construction if the API allows it
```

### Vehicle-Specific Constraint Functions

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

## API Reference

### Core Methods

#### `compute()`
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

#### Evaluation Methods
```cpp
real evalV(real s) const;      // Velocity at arc length s [m/s]
real evalAx(real s) const;     // Longitudinal acceleration at s [m/s²]
real evalAy(real s) const;     // Lateral acceleration at s [m/s²]
real evalVmax(real s) const;   // Maximum velocity at s [m/s]
```

## Tips and Best Practices

1. **Mesh Resolution**: Use appropriate spacing in your `SS_vec`. Too coarse and you'll miss important features; too fine and computation becomes slow.

2. **Initial Velocity**: Choose a reasonable initial velocity. Too high and the algorithm may fail to find a solution.

3. **Constraint Functions**: Make sure your upper and lower bounds are consistent (upper ≥ lower for all valid inputs).

4. **Debugging**: Use `get_dump()` method to identify problematic segments if the algorithm fails to converge.

5. **Performance**: The algorithm complexity scales with the number of segments, so balance accuracy with computational cost.

## Common Issues and Solutions

### Algorithm Doesn't Converge

- Check that constraint functions are well-defined for all input ranges
- Verify that the track data (curvature) is reasonable
- Reduce tolerance or increase maximum iterations of the solver
- Check for discontinuities in constraint functions

### Unrealistic Results

- Verify units (make sure everything is in SI: meters, seconds, m/s², etc.)
- Check constraint function implementations
- Validate input track data for physical reasonableness

### Performance Issues

- Reduce mesh density if acceptable for your application
- Optimize constraint function calculations
- Consider using simpler constraint models for initial testing
