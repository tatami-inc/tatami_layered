include(FetchContent)
FetchContent_Declare(
  tatami_test
  GIT_REPOSITORY https://github.com/tatami-inc/tatami_test
  GIT_TAG master
)
FetchContent_MakeAvailable(tatami_test)

find_package(ZLIB)
include(GoogleTest)

include(CheckIncludeFiles)
check_include_files(filesystem HAVE_CXX_FS)

option(DO_CODE_COVERAGE "Enable coverage testing" OFF)
if(CODE_COVERAGE AND CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    option(DO_CODE_COVERAGE ON)
endif()

# Avoid duplicating the target definition depending on whether OpenMP is to be used.
macro(create_libtest target)
  add_executable(
      ${target}
      src/convert_to_layered_sparse.cpp
      src/read_layered_sparse_from_matrix_market.cpp
      src/utils.cpp
  )

  target_link_libraries(
      ${target}
      tatami_layered
      ZLIB::ZLIB
      tatami_test
  )

  target_compile_options(${target} PRIVATE -Wall -Wextra -Wpedantic -Werror)

  if(NOT HAVE_CXX_FS)
      target_link_libraries(${target} stdc++fs)
  endif()

  if(DO_CODE_COVERAGE)
      target_compile_options(${target} PRIVATE -O0 -g --coverage)
      target_link_options(${target} PRIVATE --coverage)
  endif()

  gtest_discover_tests(${target})
endmacro()

create_libtest(libtest)

find_package(OpenMP)
if(OpenMP_FOUND)
    create_libtest(omptest)
    target_link_libraries(omptest OpenMP::OpenMP_CXX)
endif()
