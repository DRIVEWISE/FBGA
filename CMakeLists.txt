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

install(
  EXPORT FBGATargets
  FILE FBGATargets.cmake
  NAMESPACE fbga::
  DESTINATION lib/cmake/FBGA
)
