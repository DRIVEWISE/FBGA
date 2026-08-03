cmake_minimum_required(VERSION 3.24)

project(
  fbga
  LANGUAGES CXX
  VERSION 0.1.0
  HOMEPAGE_URL "https://github.com/DriveWise/FBGA"
)

include(ProjectOptions.cmake)
include(Dependencies.cmake)

add_subdirectory(src)

if(FBGA_BUILD_TESTS)
  enable_testing()
  add_subdirectory(test)
endif()

if(FBGA_BUILD_EXAMPLES)
  add_subdirectory(examples)
endif()

include(CMakePackageConfigHelpers)

configure_package_config_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/FBGAConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/FBGAConfig.cmake
  INSTALL_DESTINATION lib/cmake/FBGA
)

write_basic_package_version_file(
  ${CMAKE_CURRENT_BINARY_DIR}/FBGAConfigVersion.cmake
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion
)

install(
  EXPORT FBGATargets
  FILE FBGATargets.cmake
  NAMESPACE fbga::
  DESTINATION lib/cmake/FBGA
)

install(
  FILES
    ${CMAKE_CURRENT_BINARY_DIR}/FBGAConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/FBGAConfigVersion.cmake
  DESTINATION lib/cmake/FBGA
)
