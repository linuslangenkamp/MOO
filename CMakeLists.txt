cmake_minimum_required(VERSION 3.10)
project(moo)

set(CMAKE_CXX_STANDARD 17)

include(ExternalProject)
set(THIRD_PARTY_INSTALL_PREFIX ${CMAKE_BINARY_DIR}/third-party-install/usr/local)
set(OPENMODELICA_ROOT /home/linus/Projects/OpenModelica) # FIXME

find_package(PkgConfig)

pkg_check_modules(COINHSL coinhsl)

include(ProcessorCount)
ProcessorCount(N)

set(SOURCES
    src/base/fLGR.cpp
    src/base/trajectory.cpp
    src/base/mesh.cpp
    src/base/linalg.cpp

    src/nlp/instances/gdop/gdop.cpp
    src/nlp/instances/gdop/problem.cpp
    src/nlp/instances/gdop/gdop_strategies.cpp
    src/nlp/instances/gdop/gdop_orchestrator.cpp

    src/nlp/solvers/nlp_solver_settings.cpp
    src/nlp/nlp_scaling.cpp
    src/nlp/nlp.cpp

    src/nlp/solvers/ipopt/ipopt_adapter.cpp
    src/nlp/solvers/ipopt/ipopt_solver.cpp

    src/interfaces/openmodelica/main_opt.cpp
    src/interfaces/openmodelica/info_gdop.cpp
    src/interfaces/openmodelica/gdop_problem.cpp
    src/interfaces/openmodelica/evaluations.cpp
    src/interfaces/openmodelica/debug_om.cpp
    src/interfaces/openmodelica/sim_runtime_ext.cpp
    src/interfaces/openmodelica/strategies.cpp

    # src/nlp/instances/gdop/test_problem_impl.cpp
    # stashed/steady_state/parameter_test.cpp
)

add_library(moo SHARED ${SOURCES})

target_include_directories(moo
  PUBLIC ${CMAKE_SOURCE_DIR}/include
  PRIVATE ${CMAKE_SOURCE_DIR}/src
)

###### START FMT

option(WITH_FMT "Build and link with 'fmt' library" TRUE)

if(WITH_FMT)
  message(STATUS "Library 'fmt' will be built from source and linked (WITH_FMT=TRUE)")
  target_compile_definitions(moo PRIVATE WITH_FMT)

  # Check if the library has already been built and installed
  if(NOT EXISTS ${THIRD_PARTY_INSTALL_PREFIX}/lib/libfmt.a)
    message(STATUS "fmt library not found. Defining ExternalProject to build it.")
    ExternalProject_Add(fmt-ext
      GIT_REPOSITORY https://github.com/fmtlib/fmt.git
      GIT_TAG 11.2.0
      CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${THIRD_PARTY_INSTALL_PREFIX}
                  -DFMT_DOC=OFF
                  -DFMT_TEST=OFF
                  -DFMT_FUZZ=OFF
                  -DFMT_INSTALL=ON
                  -DCMAKE_POSITION_INDEPENDENT_CODE=ON
                  -DBUILD_SHARED_LIBS=OFF
      EXCLUDE_FROM_ALL TRUE
      BUILD_BYPRODUCTS ${THIRD_PARTY_INSTALL_PREFIX}/lib/libfmt.a
    )
    #add_dependencies(moo fmt-ext)
  else()
    message(STATUS "fmt library found. Skipping ExternalProject build.")
  endif()

  target_include_directories(moo PRIVATE
    ${THIRD_PARTY_INSTALL_PREFIX}/include
  )
  target_link_directories(moo PRIVATE
    ${THIRD_PARTY_INSTALL_PREFIX}/lib
  )
  target_link_libraries(moo PRIVATE fmt)
endif()


###### END FMT

###### START DISABLE_LOGGING

option(DISABLE_LOGGING "Disable logging macros" FALSE)

if(DISABLE_LOGGING)
  message(STATUS "Define DISABLE_LOGGING set to TRUE")
  target_compile_definitions(moo PRIVATE DISABLE_LOGGING)
else()
  message(STATUS "Define DISABLE_LOGGING set to FALSE")
endif()

###### END DISABLE_LOGGING

###### START MUMPS + Ipopt

# TODO: enable build w/o mumps or ipopt

option(WITH_MUMPS "Build and link with 'MUMPS' library" TRUE)

if(WITH_MUMPS)
  message(STATUS "Library 'MUMPS' will be built from source and linked (WITH_MUMPS=TRUE)")
  target_compile_definitions(moo PRIVATE WITH_MUMPS)
else()
  message(STATUS "Library 'MUMPS' is disabled (WITH_MUMPS=FALSE)")
endif()

# CONFIGURE_COMMAND ./configure CFLAGS=-Ofast FCFLAGS=-Ofast\ -fopenmp --enable-shared=no --enable-static=yes --prefix=${CMAKE_BINARY_DIR}/third-party-install/usr/local
ExternalProject_Add(mumps-ext
    GIT_REPOSITORY https://github.com/coin-or-tools/ThirdParty-Mumps.git
    GIT_TAG releases/3.0.10
    PATCH_COMMAND ${CMAKE_BINARY_DIR}/mumps-ext-prefix/src/mumps-ext/get.Mumps
    CONFIGURE_COMMAND ./configure --enable-shared=no --enable-static=yes --prefix=${CMAKE_BINARY_DIR}/third-party-install/usr/local
    BUILD_COMMAND make -j${N}
    BUILD_IN_SOURCE true
    EXCLUDE_FROM_ALL true
    BUILD_BYPRODUCTS ${THIRD_PARTY_INSTALL_PREFIX}/lib/libdmumps.a
)

option(WITH_IPOPT "Build and link with 'Ipopt' library" TRUE)

if(WITH_IPOPT)
  message(STATUS "Library 'Ipopt' will be built from source and linked (WITH_IPOPT=TRUE)")
  target_compile_definitions(moo PRIVATE WITH_IPOPT)
else()
  message(STATUS "Library 'Ipopt' is disabled (WITH_IPOPT=FALSE)")
endif()

ExternalProject_Add(ipopt-ext
    GIT_REPOSITORY https://github.com/coin-or/Ipopt.git
    GIT_TAG releases/3.14.19
    CONFIGURE_COMMAND ADD_CFLAGS=-Ofast ADD_CXXFLAGS=-Ofast ./configure --without-pardiso --enable-debug --enable-shared=no --enable-static=yes --prefix=${CMAKE_BINARY_DIR}/third-party-install/usr/local
    BUILD_COMMAND make -j${N}
    BUILD_IN_SOURCE true
    DEPENDS mumps-ext
    EXCLUDE_FROM_ALL true
    BUILD_BYPRODUCTS ${THIRD_PARTY_INSTALL_PREFIX}/lib/libipopt.a
)

#add_dependencies(moo ipopt-ext mumps-ext)

##### END MUMPS + Ipopt

target_include_directories(moo PUBLIC include)

target_include_directories(moo PRIVATE ${THIRD_PARTY_INSTALL_PREFIX}/include/coin-or)

target_link_directories(moo PRIVATE ${THIRD_PARTY_INSTALL_PREFIX}/lib)

#add_compile_options(-fsanitize=undefined -fsanitize=address -fstack-protector-strong)
#add_link_options(-fsanitize=undefined -fsanitize=address -fstack-protector-strong)

# TODO: Make it a clean build with -Wextra
if ("${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
    message(STATUS "Configuration: Debug")
    target_compile_options(moo PRIVATE -Wall -fno-rtti -fno-var-tracking-assignments -Wnon-virtual-dtor -pedantic -O0 -g)
    message(STATUS "Compile Options: -Wall -fno-rtti -fno-var-tracking-assignments -Wnon-virtual-dtor -pedantic -O0 -g")
    # -Wold-style-cast
else()
    message(STATUS "Configuration: Release")
    target_compile_options(moo PRIVATE -Wall -fno-rtti -fno-var-tracking-assignments -Wnon-virtual-dtor -pedantic -O3 -march=native -funroll-loops -g)
    message(STATUS "Compile Options: -Wall -fno-rtti -fno-var-tracking-assignments -Wnon-virtual-dtor -pedantic -O3 -march=native -funroll-loops -g")
endif()

target_link_libraries(moo PRIVATE ipopt coinmumps metis)
target_link_libraries(moo PRIVATE lapack blas dl gfortran gomp)
if (COINHSL_FOUND)
    message(STATUS "Coin-HSL found")
    target_link_libraries(moo PRIVATE coinhsl)
endif ()

add_subdirectory(src/modelica/example)


##### INCLUDE OM SIMULATION RUNTIME

target_link_directories(moo PRIVATE
  /home/linus/Projects/OpenModelica/build_cmake/install_cmake/lib/x86_64-linux-gnu/omc
)

target_include_directories(moo PRIVATE
  /home/linus/Projects/OpenModelica/build_cmake/install_cmake/include/omc/c
  /home/linus/Projects/OpenModelica/build_cmake/install_cmake/include/omc
  ${OPENMODELICA_ROOT}/OMCompiler
  /home/linus/Projects/OpenModelica/build_cmake/install_cmake/include/omc/suitesparse
  /home/linus/Projects/OpenModelica/build_cmake/install_cmake/include/omc/sundials
  ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(moo PRIVATE 
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

target_compile_definitions(moo PRIVATE
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
