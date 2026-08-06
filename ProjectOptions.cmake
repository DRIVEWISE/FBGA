# Language / build settings
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Build options
option(FBGA_BUILD_TESTS "Build Catch2 unit tests" ON)
option(FBGA_BUILD_EXAMPLES "Build usage examples" ON)
option(FBGA_BUILD_PYTHON "Build pybind11 module" ON)
option(FBGA_BUILD_3D "Build the fbga3d module (ported from FBGA_3D)" OFF)
option(FBGA_WARNINGS_AS_ERRORS "Treat compiler warnings as errors" OFF)

# The static libs (utils/solvers/fbga/fbga3d) get linked into the fbga_py shared module;
# on Linux, ld refuses non-PIC object code in a shared object, so every target needs PIC
# when the Python module is being built (Windows/macOS don't enforce this the same way,
# which is why this only breaks the manylinux wheel build).
if(FBGA_BUILD_PYTHON)
  set(CMAKE_POSITION_INDEPENDENT_CODE ON)
endif()

# Compiler warnings, exposed as an INTERFACE target so every in-repo target can opt
# in via `target_link_libraries(<target> PRIVATE project_warnings)` without forcing
# these flags onto fetched third-party dependencies (Catch2, cxxopts, rapidcsv, ...).
add_library(project_warnings INTERFACE)

if(MSVC)
  target_compile_options(project_warnings INTERFACE /W4)
  if(FBGA_WARNINGS_AS_ERRORS)
    target_compile_options(project_warnings INTERFACE /WX)
  endif()
else()
  target_compile_options(project_warnings INTERFACE -Wall -Wextra)
  if(FBGA_WARNINGS_AS_ERRORS)
    target_compile_options(project_warnings INTERFACE -Werror)
  endif()
endif()
