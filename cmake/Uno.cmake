target_sources(moo PRIVATE
    src/nlp/solvers/uno/adapter.cpp
    src/nlp/solvers/uno/solver.cpp
)

set(MOO_IPOPT_MUMPS_DIR "${CMAKE_CURRENT_BINARY_DIR}/third-party/Ipopt.cmake/ThirdParty/_deps/mumps-src")
set(MOO_IPOPT_MUMPS_LIB_DIR "${CMAKE_CURRENT_BINARY_DIR}/third-party/Ipopt.cmake/bin")

list(PREPEND CMAKE_INCLUDE_PATH "${MOO_IPOPT_MUMPS_DIR}/include")
list(PREPEND CMAKE_LIBRARY_PATH "${MOO_IPOPT_MUMPS_LIB_DIR}")

# Uno looks for mpiseq, Ipopt.cmake builds the same sequential MPI shim as seq.
find_library(MUMPS_MPISEQ_LIBRARY NAMES seq mpiseq PATHS "${MOO_IPOPT_MUMPS_LIB_DIR}" NO_DEFAULT_PATH)
set(MUMPS_PORD_LIBRARY "" CACHE FILEPATH "MUMPS pord library" FORCE)

set(BUILD_STATIC_LIBS ON CACHE BOOL "Build Uno static library" FORCE)
set(BUILD_SHARED_LIBS OFF CACHE BOOL "Build Uno shared library" FORCE)
set(ENABLE_TESTS OFF CACHE BOOL "Build Uno tests" FORCE)

add_subdirectory(third-party/Uno)

get_target_property(MOO_UNO_DEPENDENCY_LIBRARIES uno_dependencies INTERFACE_LINK_LIBRARIES)
if(MOO_UNO_DEPENDENCY_LIBRARIES)
    list(REMOVE_ITEM MOO_UNO_DEPENDENCY_LIBRARIES "MUMPS_PORD_LIBRARY-NOTFOUND")
    set_target_properties(uno_dependencies PROPERTIES
        INTERFACE_LINK_LIBRARIES "${MOO_UNO_DEPENDENCY_LIBRARIES}"
    )
endif()

target_include_directories(uno_dependencies INTERFACE
    "${MOO_IPOPT_MUMPS_DIR}/libseq"
)

if(TARGET uno_static)
    add_dependencies(uno_static dmumps mumps_common seq)
endif()

target_link_libraries(moo PRIVATE uno_static)
target_include_directories(moo PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/third-party/Uno)
