include(FetchContent)

FetchContent_Declare(
  eminem
  GIT_REPOSITORY https://github.com/tatami-inc/eminem
  GIT_TAG master
)

FetchContent_MakeAvailable(eminem)

FetchContent_Declare(
  tatami
  GIT_REPOSITORY https://github.com/tatami-inc/tatami
  GIT_TAG master
)

FetchContent_MakeAvailable(tatami)
