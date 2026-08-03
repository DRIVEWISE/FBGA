add_executable(
  test_fbga3d
  test_gggv_moto.cc
  test_gggv_indy.cc
  test_fbga_solver.cc
)
target_link_libraries(test_fbga3d PRIVATE fbga::fbga3d Catch2::Catch2WithMain project_warnings)

include(Catch)
# GggvIndy::setup_std() reads "./data/INDY/*.npy" relative to the working directory.
catch_discover_tests(test_fbga3d WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
