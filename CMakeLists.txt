cmake_minimum_required(VERSION 3.10)
project(OPTIMIZATION)

set(CMAKE_CXX_STANDARD 17)

file(GLOB SOURCES "src/*.cpp")

add_executable(optimization ${SOURCES})
