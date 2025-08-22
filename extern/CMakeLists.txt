include(FetchContent)

FetchContent_Declare(
  eminem
  GIT_REPOSITORY https://github.com/tatami-inc/eminem
  GIT_TAG master # ^1.2.0
)

FetchContent_Declare(
  tatami
  GIT_REPOSITORY https://github.com/tatami-inc/tatami
  GIT_TAG master # ^4.0.0
)

FetchContent_Declare(
  byteme
  GIT_REPOSITORY https://github.com/LTLA/byteme
  GIT_TAG master # ^2.0.0
)

FetchContent_Declare(
  sanisizer
  GIT_REPOSITORY https://github.com/LTLA/sanisizer
  GIT_TAG master # ^0.1.2
)

FetchContent_MakeAvailable(tatami)
FetchContent_MakeAvailable(eminem)
FetchContent_MakeAvailable(byteme)
FetchContent_MakeAvailable(sanisizer)
