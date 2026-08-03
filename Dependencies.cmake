include(FetchContent)

# Catch2 - unit testing framework, needed by test/ when FBGA_BUILD_TESTS is ON.
if(FBGA_BUILD_TESTS)
  FetchContent_Declare(
    Catch2
    GIT_REPOSITORY https://github.com/catchorg/Catch2.git
    GIT_TAG v3.15.3
  )
  FetchContent_MakeAvailable(Catch2)
  list(APPEND CMAKE_MODULE_PATH ${catch2_SOURCE_DIR}/extras)
endif()

# cxxopts / rapidcsv - needed by examples/fbga (moto example) when FBGA_BUILD_EXAMPLES is ON.
if(FBGA_BUILD_EXAMPLES)
  FetchContent_Declare(
    cxxopts
    GIT_REPOSITORY https://github.com/jarro2783/cxxopts.git
    GIT_TAG v3.3.1
  )
  FetchContent_MakeAvailable(cxxopts)

  FetchContent_Declare(
    rapidcsv
    GIT_REPOSITORY https://github.com/d99kris/rapidcsv.git
    GIT_TAG v8.99
  )
  FetchContent_MakeAvailable(rapidcsv)
endif()
