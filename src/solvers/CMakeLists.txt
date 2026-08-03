add_library(solvers)
add_library(fbga::solvers ALIAS solvers)

target_link_libraries(solvers PUBLIC fbga::utils)
target_link_libraries(solvers PRIVATE $<BUILD_INTERFACE:project_warnings>)

target_sources(
  solvers
  PRIVATE
    brent_dekker.cc
  PUBLIC
    FILE_SET HEADERS
    BASE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/include
    FILES
      include/solvers/brent_dekker.hh
)

target_compile_features(solvers PUBLIC cxx_std_20)

install(
  TARGETS solvers
  EXPORT FBGATargets
  FILE_SET HEADERS DESTINATION include
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION bin
)
