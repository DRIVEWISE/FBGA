# Language / build settings
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Build options
option(FBGA_BUILD_TESTS "Build Catch2 unit tests" ON)
option(FBGA_BUILD_EXAMPLES "Build usage examples" ON)
option(FBGA_BUILD_PYTHON "Build pybind11 module" ON)
