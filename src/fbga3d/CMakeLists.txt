add_library(fbga3d)
add_library(fbga::fbga3d ALIAS fbga3d)

target_link_libraries(fbga3d PUBLIC fbga::utils fbga::solvers)
target_link_libraries(fbga3d PRIVATE $<BUILD_INTERFACE:project_warnings> $<BUILD_INTERFACE:libnpy>)

target_sources(
  fbga3d
  PRIVATE
    gggv_moto.cc
    gggv_indy.cc
    fbga_indy.cc
    fbga_moto.cc
  PUBLIC
    FILE_SET HEADERS
    BASE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/include
    FILES
      include/fbga3d/types3d.hh
      include/fbga3d/node3d.hh
      include/fbga3d/gggv_moto.hh
      include/fbga3d/gggv_indy.hh
      include/fbga3d/fbga3d_solver.hh
      include/fbga3d/fbga_indy.hh
      include/fbga3d/fbga_moto.hh
)

target_compile_features(fbga3d PUBLIC cxx_std_20)

install(
  TARGETS fbga3d
  EXPORT FBGATargets
  FILE_SET HEADERS DESTINATION include
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION bin
)
