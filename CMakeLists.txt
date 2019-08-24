cmake_minimum_required(VERSION 3.10)
project(trng4 VERSION 4.22 DESCRIPTION "pseudo random number generator library")

include(GNUInstallDirs)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

IF (NOT CMAKE_BUILD_TYPE)
    SET(CMAKE_BUILD_TYPE Release)
ENDIF ()

add_subdirectory(trng)
add_subdirectory(examples)
