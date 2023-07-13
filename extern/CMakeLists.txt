cmake_minimum_required(VERSION 3.14 FATAL_ERROR)

include(FetchContent)

FetchContent_Declare(
  tatami
  GIT_REPOSITORY https://github.com/tatami-inc/tatami
  GIT_TAG 898039be5c0f9a4ff5de6f95ba7bdfd8d3940119
)

FetchContent_MakeAvailable(tatami)

FetchContent_Declare(
  eminem
  GIT_REPOSITORY https://github.com/tatami-inc/eminem
  GIT_TAG fc6ee36615e1d5f4a597d6ff5622960a10bc9ed7
)

FetchContent_MakeAvailable(eminem)

