add_executable(fbga_example_basic basic.cc)
target_link_libraries(fbga_example_basic PRIVATE fbga::fbga project_warnings)

add_executable(fbga_example_moto moto.cc)
target_link_libraries(fbga_example_moto PRIVATE fbga::fbga cxxopts::cxxopts rapidcsv project_warnings)
