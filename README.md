FBGA
=========

FBGA - Forward Backward Generic Acceleration constraints is a library specifically developed to solve optimal velocity profile given a complex g-g-v constraint envelope the described implementation is under review in the Robotics and Automation Letters (RAL) journal.
This library will later evolve into a generic tool including multiple dimensions constraints and a more general approach to the problem.
The library is able to compute the maximum performance longitudinal maneuver given a generic g-g-v envelope and a curvature profile (i.e., the trajectory path).

The `fbga3d` module ports the 3D forward-backward solver originally developed in [FBGA_3D](https://github.com/DRIVEWISE/FBGA_3D) into this repository (generic over the vehicle's GG-diagram model via `Fbga3dSolver<GggvModel>`, with `GggvIndy`/`GggvMoto` as the two instantiations); it's opt-in via `-DFBGA_BUILD_3D=ON` -- see below.

For usage details, please refer to the [FBGA Usage Guide](FBGA_Usage_Guide.md).

Cite 2D FBGA in your publications as follows:

```{bibtex}
@article{piazza2025real,
  title={Real-Time Velocity Profile Optimization for Time-Optimal Maneuvering With Generic Acceleration Constraints},
  author={Piazza, Mattia and Piccinini, Mattia and Taddei, Sebastiano and Biral, Francesco and Bertolazzi, Enrico},
  journal={IEEE Robotics and Automation Letters},
  volume={11},
  number={2},
  pages={1674--1681},
  year={2025},
  
  publisher={IEEE}
}
```

Cite 3D FBGA in your publications as follows:

```{bibtex}
@article{piazza2026beyond,
  title={Beyond the Apex: Online Minimum-Time Velocity Planning on Three-Dimensional Paths With Variable Friction}, 
  author={Piazza, Mattia and Langmann, Alexander and Biral, Francesco and Betz, Johannes and Piccinini, Mattia},
  journal={IEEE Robotics and Automation Letters}, 
  volume={11},
  number={8},
  pages={10034-10041},
  year={2026},
  publisher={IEEE}
}
```

<!-- download the bib -->
[__Download the bib__](bibtex.bib)



## Requirements

### General

- git
- cmake
- make

### Third party libraries

Third party libraries (Catch2, cxxopts, rapidcsv, and -- when `-DFBGA_BUILD_3D=ON` -- libnpy) are fetched automatically by CMake via `FetchContent` during configuration, pinned to fixed release tags in `Dependencies.cmake`. You need an internet connection the first time you configure the project; nothing needs to be installed manually.

## Building

```{shell}
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

For a debug build, pass `-DCMAKE_BUILD_TYPE=Debug` instead. Useful options (see `ProjectOptions.cmake`):

- `-DFBGA_BUILD_TESTS=OFF` to skip the Catch2 unit tests (`test/`).
- `-DFBGA_BUILD_EXAMPLES=OFF` to skip the usage examples (`examples/`).
- `-DFBGA_BUILD_3D=ON` to also build the `fbga3d` module (`src/fbga3d/`, `test/fbga3d/`, `examples/fbga3d/`) -- off by default.
- `-DFBGA_WARNINGS_AS_ERRORS=ON` to treat compiler warnings as errors.

### Running the tests

```{shell}
ctest --test-dir build --output-on-failure
```

### Running the examples

Examples must be run from the repository root, since they read data files with paths relative to it (e.g. `data/paper/...`, or `data/INDY/*.npy` for `fbga3d`'s `GggvIndy`):

```{shell}
./build/examples/fbga/fbga_example_basic
./build/examples/fbga/fbga_example_moto
./build/examples/fbga3d/fbga3d_example_gggv_demo         # -DFBGA_BUILD_3D=ON only
./build/examples/fbga3d/fbga3d_example_yas_marina_indy   # -DFBGA_BUILD_3D=ON only
```

## Python bindings

Python bindings (`fbga_py`, covering both `Fb2d` and `fbga3d`'s `FbgaIndy`/`FbgaMoto`) are built independently of the CMake build above, via `setuptools`/`pybind11` -- see [python_bindings/README.md](python_bindings/README.md).
