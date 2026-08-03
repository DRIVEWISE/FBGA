add_executable(fbga3d_example_gggv_demo gggv_demo.cc)
target_link_libraries(fbga3d_example_gggv_demo PRIVATE fbga::fbga3d project_warnings)

# Runs FbgaIndy against real circuit data, for comparison against FBGA_3D's own
# FBGA3D_test_INDY.cc on the same inputs -- see FBGA3D_INTEGRATION_PLAN.md.
if(TARGET rapidcsv AND TARGET cxxopts::cxxopts)
  add_executable(fbga3d_example_yas_marina_indy yas_marina_indy.cc)
  target_link_libraries(fbga3d_example_yas_marina_indy PRIVATE fbga::fbga3d rapidcsv cxxopts::cxxopts libnpy project_warnings)
endif()
