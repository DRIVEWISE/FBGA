add_library(fbga)
add_library(fbga::fbga ALIAS fbga)

target_link_libraries(fbga PUBLIC fbga::utils fbga::solvers)
target_link_libraries(fbga PRIVATE $<BUILD_INTERFACE:project_warnings>)

target_sources(
  fbga
  PRIVATE
    fb2d.cc
  PUBLIC
    FILE_SET HEADERS
    BASE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/include
    FILES
      include/fbga/fb2d.hh
)

target_compile_features(fbga PUBLIC cxx_std_20)

install(
  TARGETS fbga
  EXPORT FBGATargets
  FILE_SET HEADERS DESTINATION include
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION bin
)
