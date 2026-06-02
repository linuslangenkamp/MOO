target_sources(moo PRIVATE
    src/nlp/solvers/uno/adapter.cpp
    src/nlp/solvers/uno/solver.cpp
)

set(MOO_IPOPT_DIR "${CMAKE_CURRENT_BINARY_DIR}/third-party/Ipopt.cmake")
set(MOO_MUMPS_SRC "${MOO_IPOPT_DIR}/ThirdParty/_deps/mumps-src")

set(MUMPS_LIBRARY        dmumps       CACHE STRING "MUMPS dmumps target" FORCE)
set(MUMPS_COMMON_LIBRARY mumps_common CACHE STRING "MUMPS common target" FORCE)
set(MUMPS_MPISEQ_LIBRARY seq          CACHE STRING "MUMPS sequential MPI target" FORCE)
set(MUMPS_PORD_LIBRARY   ""           CACHE STRING "MUMPS pord target" FORCE)
set(MUMPS_INCLUDE_DIR "${MOO_MUMPS_SRC}/include" CACHE PATH "MUMPS include directory" FORCE)

set(BUILD_STATIC_LIBS ON  CACHE BOOL "Build Uno static library" FORCE)
set(BUILD_SHARED_LIBS OFF CACHE BOOL "Build Uno shared library" FORCE)
set(ENABLE_TESTS OFF      CACHE BOOL "Build Uno tests" FORCE)

add_subdirectory(third-party/Uno)

target_include_directories(uno_dependencies INTERFACE
    "${MOO_MUMPS_SRC}/libseq"
)

add_dependencies(uno_static
    dmumps
    mumps_common
    seq
)

target_link_libraries(moo PRIVATE uno_static)

target_include_directories(moo PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/third-party/Uno"
)
