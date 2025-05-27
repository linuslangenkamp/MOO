cmake_minimum_required(VERSION 3.10)
project(OPTIMIZATION)

set(CMAKE_CXX_STANDARD 17)

include(ExternalProject)
set(THIRD_PARTY_INSTALL_PREFIX ${CMAKE_BINARY_DIR}/third-party-install/usr/local)
set(OPENMODELICA_ROOT /home/linus/Projects/OpenModelica) # FIXME

find_package(PkgConfig)

pkg_check_modules(COINHSL coinhsl)

include(ProcessorCount)
ProcessorCount(N)

# CONFIGURE_COMMAND ./configure CFLAGS=-Ofast FCFLAGS=-Ofast\ -fopenmp --enable-shared=no --enable-static=yes --prefix=${CMAKE_BINARY_DIR}/third-party-install/usr/local

#[[
ExternalProject_Add(mumps-ext
    GIT_REPOSITORY https://github.com/coin-or-tools/ThirdParty-Mumps.git
    GIT_TAG releases/3.0.7
    PATCH_COMMAND ${CMAKE_BINARY_DIR}/mumps-ext-prefix/src/mumps-ext/get.Mumps
    CONFIGURE_COMMAND ./configure --enable-shared=no --enable-static=yes --prefix=${CMAKE_BINARY_DIR}/third-party-install/usr/local
    BUILD_COMMAND make -j${N}
    BUILD_IN_SOURCE true
)


ExternalProject_Add(ipopt-ext
    GIT_REPOSITORY https://github.com/coin-or/Ipopt.git
    GIT_TAG releases/3.14.16
    CONFIGURE_COMMAND ADD_CFLAGS=-Ofast ADD_CXXFLAGS=-Ofast ./configure --without-pardiso --enable-shared=no --enable-static=yes --prefix=${CMAKE_BINARY_DIR}/third-party-install/usr/local
    BUILD_COMMAND make -j${N}
    BUILD_IN_SOURCE true
    DEPENDS mumps-ext
)

add_dependencies(optimization ipopt-ext mumps-ext)

#]]

set(SOURCES
    src/base/collocation.cpp
    src/base/mesh.cpp
    src/base/linalg.cpp

    src/nlp/instances/gdop/gdop.cpp

    src/nlp/solvers/nlp_solver_flags.cpp

    src/nlp/solvers/ipopt/ipopt_adapter.cpp
    src/nlp/solvers/ipopt/ipopt_solver.cpp

    src/interfaces/openmodelica/main_opt.cpp
    src/interfaces/openmodelica/info_gdop.cpp
    src/interfaces/openmodelica/gdop_problem.cpp
    src/interfaces/openmodelica/evaluations.cpp
    src/interfaces/openmodelica/debug_om.cpp
    src/interfaces/openmodelica/sim_runtime_ext.cpp

    # src/nlp/instances/gdop/test_problem_impl.cpp
    # stashed/steady_state/parameter_test.cpp
)


add_library(optimization SHARED ${SOURCES})

target_include_directories(optimization
  PUBLIC ${CMAKE_SOURCE_DIR}/include
  PRIVATE ${CMAKE_SOURCE_DIR}/src
)


target_include_directories(optimization PUBLIC include)

target_include_directories(optimization PRIVATE ${THIRD_PARTY_INSTALL_PREFIX}/include/coin-or)

target_link_directories(optimization PRIVATE ${THIRD_PARTY_INSTALL_PREFIX}/lib)

#add_compile_options(-fsanitize=undefined -fsanitize=address -fstack-protector-strong)
#add_link_options(-fsanitize=undefined -fsanitize=address -fstack-protector-strong)

if (${CMAKE_BUILD_TYPE} STREQUAL "Debug")
    message(STATUS "Configuration: Debug")
    target_compile_options(optimization PRIVATE -Wall -fno-rtti -fno-var-tracking-assignments -O0 -g)
else()
    message(STATUS "Configuration: Release")
    target_compile_options(optimization PRIVATE -Wall -fno-rtti -fno-var-tracking-assignments -O2 -g)
endif()

target_link_libraries(optimization PRIVATE ipopt coinmumps metis)
target_link_libraries(optimization PRIVATE lapack blas dl gfortran gomp)
if (COINHSL_FOUND)
    message(STATUS "Coin-HSL found")
    target_link_libraries(optimization PRIVATE coinhsl)
endif ()

add_subdirectory(src/modelica/example)


##### INCLUDE OM SIMULATION RUNTIME

target_link_directories(optimization PRIVATE
  /home/linus/Projects/OpenModelica/build_cmake/install_cmake/lib/x86_64-linux-gnu/omc
)

target_include_directories(optimization PRIVATE
  /home/linus/Projects/OpenModelica/build_cmake/install_cmake/include/omc/c
  /home/linus/Projects/OpenModelica/build_cmake/install_cmake/include/omc
  ${OPENMODELICA_ROOT}/OMCompiler
  /home/linus/Projects/OpenModelica/build_cmake/install_cmake/include/omc/suitesparse
  /home/linus/Projects/OpenModelica/build_cmake/install_cmake/include/omc/sundials
  ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(optimization PRIVATE 
    SimulationRuntimeC
    OpenModelicaRuntimeC
    omcgc
    zlib
    lapack
    blas
    gfortran
    pthread
    dl
    m
    stdc++
)

target_compile_definitions(optimization PRIVATE
  OM_HAVE_PTHREADS
  OMC_DLL_MAIN_DEFINE
  OPENMODELICA_XML_FROM_FILE_AT_RUNTIME
  OMC_MODEL_PREFIX=include_test
  OMC_NUM_MIXED_SYSTEMS=0
  OMC_NUM_LINEAR_SYSTEMS=0
  OMC_NUM_NONLINEAR_SYSTEMS=0
  OMC_NDELAY_EXPRESSIONS=0
  OMC_NVAR_STRING=0
)

##### INCLUDE OM SIMULATION RUNTIME
