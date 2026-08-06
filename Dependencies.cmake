include(FetchContent)

# Catch2 - unit testing framework, needed by test/ when FBGA_BUILD_TESTS is ON.
if(FBGA_BUILD_TESTS)
  FetchContent_Declare(
    Catch2
    GIT_REPOSITORY https://github.com/catchorg/Catch2.git
    GIT_TAG v3.15.3
  )
  FetchContent_MakeAvailable(Catch2)
  list(APPEND CMAKE_MODULE_PATH ${catch2_SOURCE_DIR}/extras)
endif()

# cxxopts / rapidcsv - needed by examples/fbga (moto example) when FBGA_BUILD_EXAMPLES is ON.
if(FBGA_BUILD_EXAMPLES)
  FetchContent_Declare(
    cxxopts
    GIT_REPOSITORY https://github.com/jarro2783/cxxopts.git
    GIT_TAG v3.3.1
  )
  FetchContent_MakeAvailable(cxxopts)

  FetchContent_Declare(
    rapidcsv
    GIT_REPOSITORY https://github.com/d99kris/rapidcsv.git
    GIT_TAG v8.99
  )
  FetchContent_MakeAvailable(rapidcsv)
endif()

# libnpy - header-only .npy reader, needed by fbga3d's GggvIndy (loads spline data
# from .npy files) when FBGA_BUILD_3D is ON. Upstream ships no CMakeLists.txt (meson
# only), so FetchContent_MakeAvailable() just populates the source; wrap the header in
# an INTERFACE target ourselves.
if(FBGA_BUILD_3D OR FBGA_BUILD_PYTHON)
  FetchContent_Declare(
    libnpy
    GIT_REPOSITORY https://github.com/llohse/libnpy.git
    GIT_TAG v1.0.1
  )
  FetchContent_MakeAvailable(libnpy)

  add_library(libnpy INTERFACE)
  target_include_directories(libnpy INTERFACE ${libnpy_SOURCE_DIR}/include)
endif()

# pybind11 - needed by src/python when FBGA_BUILD_PYTHON is ON. find_package first
# (scikit-build-core sets pybind11_DIR via its pybind11 build dependency), FetchContent
# fallback for a plain CMake configure.
if(FBGA_BUILD_PYTHON)
  find_package(pybind11 CONFIG QUIET)
  if(NOT pybind11_FOUND)
    FetchContent_Declare(
      pybind11
      GIT_REPOSITORY https://github.com/pybind/pybind11.git
      GIT_TAG v2.13.6
    )
    FetchContent_MakeAvailable(pybind11)
  endif()
endif()
