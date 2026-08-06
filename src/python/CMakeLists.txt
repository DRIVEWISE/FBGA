pybind11_add_module(fbga_py bindings.cc)

target_link_libraries(fbga_py PRIVATE fbga::utils fbga::solvers fbga::fbga fbga::fbga3d)

target_compile_features(fbga_py PRIVATE cxx_std_20)

install(TARGETS fbga_py LIBRARY DESTINATION . COMPONENT python_module)
