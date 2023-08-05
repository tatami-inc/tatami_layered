include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/afd902e992b720d1b3e106bc5e425a5768872265.zip
)

# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

# Avoid installing GoogleTest when installing this project.
option(INSTALL_GTEST "Enable installation of googletest." OFF)

FetchContent_MakeAvailable(googletest)

enable_testing()

find_package(ZLIB)

include(GoogleTest)

include(CheckIncludeFiles)
check_include_files(filesystem HAVE_CXX_FS)

set(CODE_COVERAGE OFF CACHE BOOL "Enable coverage testing")

# Avoid duplicating the target definition depending on whether OpenMP is to be used.
macro(create_libtest target)
  add_executable(
      ${target}
      src/convert_to_layered_sparse.cpp
      src/read_layered_sparse_from_matrix_market.cpp
  )
  
  target_link_libraries(
      ${target}
      gtest_main
      tatami_layered
      ZLIB::ZLIB
  )

  if(NOT HAVE_CXX_FS)
      target_link_libraries(${target} stdc++fs)
  endif()

  if(CODE_COVERAGE AND CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
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
