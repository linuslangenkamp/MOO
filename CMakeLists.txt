cmake_minimum_required(VERSION 3.10)
project(OPTIMIZATION)

set(CMAKE_CXX_STANDARD 17)

set(SOURCES
    src/collocation.cpp
    src/main.cpp
    src/mesh.cpp
    src/nlp.cpp
    src/problem.cpp
    src/base/linalg.cpp
)

include_directories(${CMAKE_SOURCE_DIR}/src ${CMAKE_SOURCE_DIR}/src/base)

add_executable(optimization ${SOURCES})