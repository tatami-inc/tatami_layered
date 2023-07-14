cmake_minimum_required(VERSION 3.14 FATAL_ERROR)

include(FetchContent)

FetchContent_Declare(
  eminem
  GIT_REPOSITORY https://github.com/tatami-inc/eminem
  GIT_TAG 2606bcd44b8901c869deaf354f4b1df5a280eb89
)

FetchContent_MakeAvailable(eminem)

FetchContent_Declare(
  tatami
  GIT_REPOSITORY https://github.com/tatami-inc/tatami
  GIT_TAG 898039be5c0f9a4ff5de6f95ba7bdfd8d3940119
)

FetchContent_MakeAvailable(tatami)
