add_executable(test_fbga test_fb2d.cc)
target_link_libraries(test_fbga PRIVATE fbga::fbga Catch2::Catch2WithMain project_warnings)

include(Catch)
catch_discover_tests(test_fbga)
