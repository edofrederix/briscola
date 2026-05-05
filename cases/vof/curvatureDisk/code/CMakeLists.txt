cmake_minimum_required(VERSION 3.22)
project(code CXX)

if(NOT DEFINED ENV{BRISCOLA} OR "$ENV{BRISCOLA}" STREQUAL "")
    message(FATAL_ERROR "BRISCOLA environment variable not set")
endif()

include($ENV{BRISCOLA}/cmake/Briscola.cmake)

find_package(Briscola REQUIRED)

add_library(code SHARED curvatureDisk.C)
target_link_libraries(code PRIVATE Briscola)
