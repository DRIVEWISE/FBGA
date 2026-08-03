add_library(utils)
add_library(fbga::utils ALIAS utils)

target_sources(
  utils
  PRIVATE
    gg_utils.cc
  PUBLIC
    FILE_SET HEADERS
    BASE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/include
    FILES
      include/utils/types.hh
      include/utils/gg_utils.hh
      include/utils/segment.hh
)

target_compile_features(utils PUBLIC cxx_std_20)

install(
  TARGETS utils
  EXPORT FBGATargets
  FILE_SET HEADERS DESTINATION include
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION bin
)
