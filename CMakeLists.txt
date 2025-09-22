cmake_minimum_required(VERSION 3.10)
project(moo LANGUAGES C CXX Fortran)

set(CMAKE_CXX_STANDARD 17)

### Options

set(MOO_MAKEPROGRAM "" CACHE STRING "Make program used to build third party dependencies of MOO library.")
set(MOO_SUBMAKE_FLAGS "" CACHE STRING "Extra flags passed to sub-make when building MUMPS/Ipopt (leave empty for standard hierarchy, set e.g. -j{PROCS} for parallel build, or -j1 for single-threaded). This has no effect if CMake generator is Unix or MSYS Makefiles!")

set(MOO_LAPACK_LIB "" CACHE STRING "Override path to LAPACK library")
set(MOO_GFORTRAN_LIB "" CACHE STRING "Override path to gfortran library")

option(MOO_WITH_RADAU "Enable RADAU build and linking" ON)
option(MOO_WITH_GDOPT "Build MOO with GDOPT support" ON)
option(MOO_WITH_GDOPT_EXAMPLE "Build MOO with GDOPT example executable" OFF)

##### START setup

message(STATUS "C Compiler: ${CMAKE_C_COMPILER}")
message(STATUS "C++ Compiler: ${CMAKE_CXX_COMPILER}")
message(STATUS "Fortran Compiler: ${CMAKE_Fortran_COMPILER}")

if(WIN32)
    set(LIB_EXT .dll)
elseif(APPLE)
    set(LIB_EXT .dylib)
else()
    set(LIB_EXT .so)
endif()

message(STATUS "Target library extension: ${LIB_EXT}")

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    set(IS_GCC TRUE)
else()
    set(IS_GCC FALSE)
endif()

##### END setup

##### END make

include(ExternalProject)
set(THIRD_PARTY_INSTALL_PREFIX ${CMAKE_CURRENT_BINARY_DIR}/third-party-install/usr/local)

find_package(PkgConfig)

###### START DEPENDENCIES

pkg_check_modules(COINHSL coinhsl)

set(MOO_EXTRA_LIBS "")

# LAPACK
if(MOO_LAPACK_LIB)
    message(STATUS "Using user-provided LAPACK: ${MOO_LAPACK_LIB}")
    list(APPEND MOO_EXTRA_LIBS ${MOO_LAPACK_LIB})
else()
    message(STATUS "MOO_LAPACK_LIB not defined: searching for LAPACK")
    find_package(LAPACK REQUIRED)

    if(LAPACK_FOUND)
        set(MOO_LAPACK_LIB ${LAPACK_LIBRARIES})
        message(STATUS "Found LAPACK: ${MOO_LAPACK_LIB}")
        list(APPEND MOO_EXTRA_LIBS ${MOO_LAPACK_LIB})
    else()
        message(FATAL_ERROR "LAPACK library not found")
    endif()
endif()

string(JOIN " " MOO_LAPACK_LIB_STRING ${MOO_LAPACK_LIB})

# gfortran runtime
if(MOO_GFORTRAN_LIB)
    message(STATUS "Using user-provided gfortran: ${MOO_GFORTRAN_LIB}")
    list(APPEND MOO_EXTRA_LIBS ${MOO_GFORTRAN_LIB})
else()
    message(STATUS "MOO_GFORTRAN_LIB not defined: searching for gfortran")
    find_library(MOO_GFORTRAN_LIB NAMES gfortran)

    if(NOT MOO_GFORTRAN_LIB)
        message(STATUS "gfortran library not found in default paths. Attempting to find via the compiler.")
        find_program(GFORTRAN_EXECUTABLE NAMES gfortran)

        if(GFORTRAN_EXECUTABLE)
            message(STATUS "Found gfortran executable: ${GFORTRAN_EXECUTABLE}")

            execute_process(
                COMMAND "${GFORTRAN_EXECUTABLE}" -print-file-name=libgfortran${LIB_EXT}
                OUTPUT_VARIABLE GFORTRAN_LIB_PATH
                OUTPUT_STRIP_TRAILING_WHITESPACE
                RESULT_VARIABLE EXECUTE_RESULT
            )

            if(GFORTRAN_LIB_PATH)
                message(STATUS "Found gfortran library path using the compiler: ${GFORTRAN_LIB_PATH}")
                set(MOO_GFORTRAN_LIB "${GFORTRAN_LIB_PATH}")
            else()
                message(STATUS "Failed to get gfortran library path from the compiler.")
            endif()
        else()
            message(STATUS "gfortran executable not found.")
        endif()
    else()
        message(STATUS "Found gfortran: ${MOO_GFORTRAN_LIB}")
        list(APPEND MOO_EXTRA_LIBS ${MOO_GFORTRAN_LIB})
    endif()

    if(NOT MOO_GFORTRAN_LIB)
        message(WARNING "gfortran library not found: linking may still work if system finds it automatically")
    endif()
endif()

###### END DEPENDENCIES

###### START BUILD + COMPILER FLAGS

# Attention:
# force release build if included in 3rdParty tool - explicitly set MOO_BUILD_TYPE to debug for debug build!
if(NOT DEFINED MOO_BUILD_TYPE)
  set(MOO_BUILD_TYPE "RelWithDebInfo")
  message(STATUS "MOO_BUILD_TYPE not defined: defaulting to RelWithDebInfo")
endif()

if(MOO_NATIVE_BUILD)
    set(MOO_O3_BUILD_FLAGS
        -O3
        -march=native
    )
    message(STATUS "Adding native CPU optimizations: ${MOO_O3_BUILD_FLAGS}")

else()
    set(MOO_O3_BUILD_FLAGS -O3)
    message(STATUS "MOO_NATIVE_BUILD not defined: adding no additional compiler flags")
endif()

set(MOO_MUMPS_FORTRAN_FLAGS -O3)

string(JOIN " " MOO_MUMPS_FORTRAN_FLAGS_STRING ${MOO_MUMPS_FORTRAN_FLAGS})

###### END BUILD + COMPILER FLAGS

add_library(moo SHARED)

target_sources(moo PRIVATE
    # base sources for GDOP and similar problems
    src/base/fLGR.cpp
    src/base/trajectory.cpp
    src/base/mesh.cpp
    src/base/log.cpp
    src/base/linalg.cpp

    # General Dynamic Optimization Problem
    src/nlp/instances/gdop/gdop.cpp
    src/nlp/instances/gdop/problem.cpp
    src/nlp/instances/gdop/strategies.cpp
    src/nlp/instances/gdop/orchestrator.cpp

    # generic NLP
    src/nlp/solvers/nlp_solver_settings.cpp
    src/nlp/nlp_scaling.cpp
    src/nlp/nlp.cpp

    # Ipopt Interface for generic NLP
    src/nlp/solvers/ipopt/adapter.cpp
    src/nlp/solvers/ipopt/solver.cpp

    # C interface
    src/interfaces/c/problem.cpp

    # Sources to call simulations for C interface
    src/simulation/integrator/integrator_util.cpp
    src/simulation/integrator/integrator.cpp
    src/simulation/radau/radau_builder.cpp
    src/simulation/radau/radau_integrator.cpp
    src/simulation/radau/test.cpp
)

if (MOO_WITH_GDOPT)
    message(STATUS "Building MOO with GDOPT support")
    target_sources(moo PRIVATE src/interfaces/gdopt/main_gdopt.cpp )
endif ()

target_include_directories(moo
  PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src
)

###### START FMT

add_subdirectory(third-party/fmt)
target_compile_options(fmt PRIVATE -fPIC)

###### END FMT

###### START MUMPS + Ipopt

add_subdirectory(third-party/Ipopt.cmake)

target_include_directories(moo PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/third-party/Ipopt.cmake/Ipopt/sources/src/Interfaces
    ${CMAKE_CURRENT_SOURCE_DIR}/third-party/Ipopt.cmake/Ipopt/sources/src/Common
    ${CMAKE_CURRENT_SOURCE_DIR}/third-party/Ipopt.cmake/Ipopt/sources/src/Algorithm
    ${CMAKE_CURRENT_SOURCE_DIR}/third-party/Ipopt.cmake/Ipopt/sources/src/LinAlg
)

###### END MUMPS + Ipopt

###### START RADAU5

if(MOO_WITH_RADAU)
    add_subdirectory(third-party/RADAU)
endif()

##### END RADAU5

target_include_directories(moo PUBLIC include)
target_include_directories(moo PUBLIC src)

target_link_directories(moo PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/Dependencies/lib)

# add_compile_options(-fsanitize=undefined -fsanitize=address -fstack-protector-strong)
# add_link_options(-fsanitize=undefined -fsanitize=address -fstack-protector-strong)

# TODO: Make it a clean build with -Wextra

set(MOO_COMPILE_OPTIONS -Wall -fno-rtti -Wnon-virtual-dtor -pedantic)

if(IS_GCC AND NOT "${MOO_BUILD_TYPE}" STREQUAL "Debug")
    list(APPEND MOO_COMPILE_OPTIONS -fno-var-tracking-assignments)
endif()

if("${MOO_BUILD_TYPE}" STREQUAL "Debug")
    message(STATUS "Configuration: Debug")
    list(APPEND MOO_COMPILE_OPTIONS -O0 -g)
else()
    message(STATUS "Configuration: RelWithDebInfo / Release")
    list(APPEND MOO_COMPILE_OPTIONS -g ${MOO_O3_BUILD_FLAGS})
endif()

target_compile_options(moo PRIVATE ${MOO_COMPILE_OPTIONS})
message(STATUS "Compile Options: ${MOO_COMPILE_OPTIONS}")

target_link_libraries(moo PUBLIC fmt)
target_link_libraries(moo PRIVATE ipopt)
target_link_libraries(moo PRIVATE ${MOO_EXTRA_LIBS})

if(MOO_WITH_RADAU)
    target_link_libraries(moo PRIVATE radau)
endif()

if (COINHSL_FOUND)
    target_link_libraries(moo PRIVATE coinhsl)
endif ()

# hide symbols
if(WIN32)
    target_compile_definitions(moo PRIVATE MOO_DLL_EXPORT)
elseif (${CMAKE_SYSTEM_NAME} STREQUAL "Darwin")
    target_link_libraries(moo PRIVATE dl)
    target_link_directories(moo PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/third-party/mumps/.libs/)
    target_link_directories(moo PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/third-party/ipopt/src/.libs/)
    target_link_options(moo PRIVATE -Wl,-hidden-lcoinmumps -Wl,-hidden-lipopt)
    target_compile_options(moo PRIVATE -fvisibility=hidden )
else()
    target_link_libraries(moo PRIVATE dl)
    target_link_options(moo PRIVATE -Wl,--exclude-libs,ALL)
    target_compile_options(moo PRIVATE -fvisibility=hidden)
endif()

install(TARGETS moo
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})

if(MOO_WITH_GDOPT_EXAMPLE)
    add_executable(gdopt src/runner/entry.cpp src/runner/generated.c)
    target_include_directories(gdopt PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/include
        ${CMAKE_CURRENT_SOURCE_DIR}/src
    )
    target_link_libraries(gdopt PRIVATE moo)
    message(STATUS "Building example executable ./gdopt")
endif()
