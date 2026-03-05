FBGA
=========

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

Some third party libraries are needed to compile the project.

You will need to run the following command to install the third party libraries:

```{shell}
./third_party.sh
```

The library is self-contained and standalone, so you don't need to install any other library to use it as is. Third party libraries are used for testing and development purposes.

<!-- Note -->

**Note**: The third party libraries are not included in the repository, so you need to run the command above to download them. Please make sure you have an internet connection when running the command and that you have the necessary permissions to install the libraries.

## Building

To install the project you need to run the following commands:

```{shell}
./build.sh
```

Default will build the project in release mode, if you want to build in debug mode you can run:

```{shell}
./build.sh -debug
```

For further information you can run:

```{shell}
./build.sh -h
```

### Alternative build

If you want to install the project in a different directory you can run:

```{shell}
mkdir -p build
cd build
cmake -DCMAKE_BUILD_TYPE="$BUILD_TYPE" ..
make -j
```

### Running the examples 

To run the tests you can use the following command:

```{shell}
./bin/GIGI_test_moto.exe
```
