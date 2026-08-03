add_executable(test_solvers test_brent_dekker.cc)
target_link_libraries(test_solvers PRIVATE fbga::solvers Catch2::Catch2WithMain)

include(Catch)
catch_discover_tests(test_solvers)
