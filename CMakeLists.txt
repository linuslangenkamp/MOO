cmake_minimum_required(VERSION 3.10)
project(OPTIMIZATION)

set(CMAKE_CXX_STANDARD 17)

include(ExternalProject)
set(THIRD_PARTY_INSTALL_PREFIX ${CMAKE_BINARY_DIR}/third-party-install/usr/local)

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

#--enable-debug 
ExternalProject_Add(ipopt-ext
    GIT_REPOSITORY https://github.com/coin-or/Ipopt.git
    GIT_TAG releases/3.14.16
    CONFIGURE_COMMAND ADD_CFLAGS=-Ofast ADD_CXXFLAGS=-Ofast ./configure --without-pardiso --enable-shared=no --enable-static=yes --prefix=${CMAKE_BINARY_DIR}/third-party-install/usr/local
    BUILD_COMMAND make -j${N}
    BUILD_IN_SOURCE true
    DEPENDS mumps-ext
)
#]]

set(SOURCES
    src/base/collocation.cpp
    src/main.cpp
    src/base/mesh.cpp
    src/impl/gdop/gdop.cpp
    src/base/linalg.cpp
    src/interfaces/ipopt_adapter.cpp
    src/interfaces/ipopt_solver.cpp
    src/impl/gdop/test_problem_impl.cpp
)

include_directories(${CMAKE_SOURCE_DIR}/src ${CMAKE_SOURCE_DIR}/src/base ${CMAKE_SOURCE_DIR}/src/interfaces)

add_executable(optimization ${SOURCES})

target_include_directories(optimization PRIVATE ${THIRD_PARTY_INSTALL_PREFIX}/include/coin-or)

target_link_directories(optimization PRIVATE ${THIRD_PARTY_INSTALL_PREFIX}/lib)

#add_dependencies(optimization ipopt-ext mumps-ext)

target_compile_options(optimization PRIVATE -Wall -Werror -fno-rtti -fno-var-tracking-assignments -Ofast -fopenmp)
target_link_libraries(optimization PRIVATE ipopt coinmumps metis)
target_link_libraries(optimization PRIVATE lapack blas dl gfortran gomp)
if (COINHSL_FOUND)
    target_link_libraries(optimization PRIVATE coinhsl)
endif ()
