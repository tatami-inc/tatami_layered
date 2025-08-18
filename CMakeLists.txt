cmake_minimum_required(VERSION 3.24)

project(tatami_layered
    VERSION 2.0.0
    DESCRIPTION "Create layered tatami matrices"
    LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)

add_library(tatami_layered INTERFACE)
add_library(tatami::tatami_layered ALIAS tatami_layered)

option(TATAMI_LAYERED_FETCH_EXTERN "Automatically fetch tatami_layered's external dependencies." ON)
if(TATAMI_LAYERED_FETCH_EXTERN)
    add_subdirectory(extern)
else()
    find_package(tatami_tatami 4.0.0 CONFIG REQUIRED)
    find_package(tatami_eminem 1.2.0 CONFIG REQUIRED)
    find_package(ltla_byteme 2.0.0 CONFIG REQUIRED)
    find_package(ltla_sanisizer 0.1.0 CONFIG REQUIRED)
endif()

target_link_libraries(tatami_layered INTERFACE tatami::tatami tatami::eminem ltla::byteme ltla::sanisizer)

# Switch between include directories depending on whether the downstream is
# using the build directly or is using the installed package.
include(GNUInstallDirs)
target_include_directories(tatami_layered
    INTERFACE 
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/tatami_tatami_layered>"
)

# Building the test-related machinery, if we are compiling this library directly.
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    option(TATAMI_LAYERED_TESTS "Build tatami_layered's test suite." ON)
else()
    option(TATAMI_LAYERED_TESTS "Build tatami_layered's test suite." OFF)
endif()

if(TATAMI_LAYERED_TESTS)
    include(CTest)
    if(BUILD_TESTING)
        add_subdirectory(tests)
    endif()
endif()

# Installing for find_package.
include(CMakePackageConfigHelpers)

install(DIRECTORY include/
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/tatami_tatami_layered)

install(TARGETS tatami_layered
    EXPORT tatami_layeredTargets)

install(EXPORT tatami_layeredTargets
    FILE tatami_tatami_layeredTargets.cmake
    NAMESPACE tatami::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/tatami_tatami_layered)

configure_package_config_file(${CMAKE_CURRENT_SOURCE_DIR}/cmake/Config.cmake.in
    "${CMAKE_CURRENT_BINARY_DIR}/tatami_tatami_layeredConfig.cmake"
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/tatami_tatami_layered)

install(FILES "${CMAKE_CURRENT_BINARY_DIR}/tatami_tatami_layeredConfig.cmake"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/tatami_tatami_layered)
