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

add_executable(
  libtest
  src/convert_to_layered_sparse.cpp
  src/read_layered_sparse_from_matrix_market.cpp
  src/utils.cpp
)

target_link_libraries(
  libtest
  tatami_layered
  ZLIB::ZLIB
  tatami_test
)

target_compile_options(libtest PRIVATE -Wall -Wextra -Wpedantic -Werror)

if(NOT HAVE_CXX_FS)
  target_link_libraries(libtest stdc++fs)
endif()

option(CODE_COVERAGE "Enable coverage testing" OFF)
if(CODE_COVERAGE AND CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  target_compile_options(libtest PRIVATE -O0 -g --coverage)
  target_link_options(libtest PRIVATE --coverage)
endif()

option(FORCE_COPY "Force copy in tatami::fetch" OFF)
if (FORCE_COPY)
  target_compile_definitions(libtest PRIVATE TATAMI_DEBUG_FORCE_COPY=1)
endif()
