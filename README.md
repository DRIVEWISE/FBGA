FBGA
=========

For FBGA 3D road refer to repository [FBGA_3D](https://github.com/DRIVEWISE/FBGA_3D). In the future, the FBGA library will be merged with the FBGA_3D repository to provide a single library for both 2D and 3D road scenarios.

FBGA - Forward Backward Generic Acceleration constraints is a library specifically developed to solve optimal velocity profile given a complex g-g-v constraint envelope the described implementation is under review in the Robotics and Automation Letters (RAL) journal.
This library will later evolve into a generic tool including multiple dimensions constraints and a more general approach to the problem.
The library is able to compute the maximum performance longitudinal maneuver given a generic g-g-v envelope and a curvature profile (i.e., the trajectory path).

For usage details, please refer to the [FBGA Usage Guide](FBGA_Usage_Guide.md).

Cite as:

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
<!-- download the bib -->
[__Download the bib__](bibtex.bib)



## Requirements

### General

- git
- cmake
- make

### Third party libraries

Third party libraries (Catch2, cxxopts, rapidcsv) are fetched automatically by CMake via `FetchContent` during configuration, pinned to fixed release tags in `Dependencies.cmake`. You need an internet connection the first time you configure the project; nothing needs to be installed manually.

## Building

```{shell}
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

For a debug build, pass `-DCMAKE_BUILD_TYPE=Debug` instead. Useful options (see the top-level `CMakeLists.txt`):

- `-DFBGA_BUILD_TESTS=OFF` to skip the Catch2 unit tests (`test/`).
- `-DFBGA_BUILD_EXAMPLES=OFF` to skip the usage examples (`examples/`).
- `-DFBGA_BUILD_PYTHON=OFF` to skip the Python bindings.

### Running the tests

```{shell}
ctest --test-dir build --output-on-failure
```

### Running the examples

Examples must be run from the repository root, since they read data files with paths relative to it (e.g. `data/paper/...`):

```{shell}
./build/examples/fbga/fbga_example_basic
./build/examples/fbga/fbga_example_moto
```
