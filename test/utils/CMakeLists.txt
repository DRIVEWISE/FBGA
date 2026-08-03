add_executable(
  test_utils
  test_gg_utils.cc
  test_segment.cc
)
target_link_libraries(test_utils PRIVATE fbga::utils Catch2::Catch2WithMain project_warnings)

include(Catch)
catch_discover_tests(test_utils)
