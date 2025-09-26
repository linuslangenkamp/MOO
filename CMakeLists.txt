cmake_minimum_required(VERSION 3.10) # 3.21 for FindMsys
project(moo LANGUAGES C CXX Fortran)

set(CMAKE_CXX_STANDARD 17)

### Options

set(MOO_MAKEPROGRAM "" CACHE STRING "Make program used to build third party dependencies of MOO library.")
set(MOO_SUBMAKE_FLAGS "" CACHE STRING "Extra flags passed to sub-make when building MUMPS/Ipopt (leave empty for standard hierarchy, set e.g. -j{PROCS} for parallel build, or -j1 for single-threaded). This has no effect if CMake generator is Unix or MSYS Makefiles!")

set(MOO_LAPACK_LIB "" CACHE STRING "Override path to LAPACK library")
set(MOO_METIS_LIB "" CACHE STRING "Override path to METIS library")
set(MOO_GFORTRAN_LIB "" CACHE STRING "Override path to gfortran library")

option(MOO_WITH_RADAU "Enable RADAU build and linking" ON)
option(MOO_WITH_GDOPT "Build MOO with GDOPT support" ON)

option(MOO_TESTS "Build MOO with test executables" OFF)

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

set(MOO_IS_MAKEFILE_GENERATOR OFF)
if(CMAKE_GENERATOR STREQUAL "Unix Makefiles" OR CMAKE_GENERATOR STREQUAL "MSYS Makefiles")
    set(MOO_IS_MAKEFILE_GENERATOR ON)
endif()

##### END setup

##### START make

# A make program is necessary for the legacy configure builds of Ipopt and MUMPS

set(MOO_MAKEPROGRAM_INTERNAL "")
set(MOO_SUBMAKE_FLAGS_INTERNAL "")
set(MOO_SUBMAKE_CMD "")

# determine the make program if not already defined
if(MOO_MAKEPROGRAM STREQUAL "")
    if(MOO_IS_MAKEFILE_GENERATOR)
        # Unix/MSYS Makefile generators support $(MAKE) environment variable -> needed for parallel jobserver makes
        set(MOO_MAKEPROGRAM_INTERNAL "\$(MAKE)")
        message(STATUS "Using automatic parallel jobserver make (generator supports $(MAKE)).")
    # fallbacks cause MOO_MAKEPROGRAM not set
    elseif(DEFINED ENV{MAKE})
        # use MAKE from environment
        set(MOO_MAKEPROGRAM_INTERNAL "$ENV{MAKE}")
        message(STATUS "Detected MAKE from environment: using ${MOO_MAKEPROGRAM_INTERNAL}")
    else()
        # find in PATH
        message(STATUS "Attempting to find a make program in PATH...")
        find_program(MAKE_PROGRAM "make")
        set(MOO_MAKEPROGRAM_INTERNAL "${MAKE_PROGRAM}")
        if(MAKE_PROGRAM STREQUAL "MAKE_PROGRAM-NOTFOUND")
            set(MOO_MAKEPROGRAM_INTERNAL "\$(MAKE)")
            message(WARNING "Could not find a make program in PATH. Falling back to $(MAKE) environment variable. Please set MOO_MAKEPROGRAM to specify an appropriate make program.")
            if(CMAKE_MAKE_PROGRAM MATCHES "ninja$")
                message(SEND_ERROR "Ninja build does not support reading $(MAKE) environment variable. Build will likely fail.")
            endif()
        else()
            message(STATUS "Found make program in PATH: ${MAKE_PROGRAM}")
        endif()
    endif()
else()
    set(MOO_MAKEPROGRAM_INTERNAL "${MOO_MAKEPROGRAM}")
    message(STATUS "Using user specified MOO_MAKEPROGRAM: ${MOO_MAKEPROGRAM_INTERNAL}.")
endif()

# determine sub-make flags only if the generator does not support automatic jobserver
if(NOT MOO_IS_MAKEFILE_GENERATOR)
    # handle parallelism
    if(MOO_SUBMAKE_FLAGS STREQUAL "")
        message(WARNING "Cannot use automatic parallel make jobserver, because generator does not support \$(MAKE). "
                        "MOO_SUBMAKE_FLAGS not set: using default resolution order (user > MAKEFLAGS > half of available cores). "
                        "For a safe single-threaded build, set MOO_SUBMAKE_FLAGS=\"\" explicitly.")
        if(DEFINED ENV{MAKEFLAGS} AND NOT "$ENV{MAKEFLAGS}" STREQUAL "")
            set(MOO_SUBMAKE_FLAGS_INTERNAL "$ENV{MAKEFLAGS}")
            message(STATUS "Detected MAKEFLAGS from environment: ${MOO_SUBMAKE_FLAGS_INTERNAL}")
        else()
            include(ProcessorCount)
            ProcessorCount(NPROC)
            math(EXPR SUBMAKE_JOBS "${NPROC}/2")
            if(SUBMAKE_JOBS LESS 1)
                set(SUBMAKE_JOBS 1)
            endif()
            set(MOO_SUBMAKE_FLAGS_INTERNAL "-j${SUBMAKE_JOBS}")
            message(STATUS "Using default sub-make parallelism (half of available cores, ${MOO_SUBMAKE_FLAGS_INTERNAL})")
        endif()
    else()
        set(MOO_SUBMAKE_FLAGS_INTERNAL "${MOO_SUBMAKE_FLAGS}")
        message(STATUS "Using user-specified sub-make flags: ${MOO_SUBMAKE_FLAGS_INTERNAL}")
    endif()
endif()

set(MOO_SUBMAKE_CMD ${MOO_MAKEPROGRAM_INTERNAL} ${MOO_SUBMAKE_FLAGS_INTERNAL})
message(STATUS "Make command used to build third-party dependencies of MOO library: ${MOO_SUBMAKE_CMD}")

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

# METIS
if(MOO_METIS_LIB)
    message(STATUS "Using user-provided METIS: ${MOO_METIS_LIB}")
    list(APPEND MOO_EXTRA_LIBS ${MOO_METIS_LIB})
else()
    message(STATUS "MOO_METIS_LIB not defined: searching for METIS")
    find_library(MOO_METIS_LIB NAMES metis)
    if(MOO_METIS_LIB)
        message(STATUS "Found METIS: ${MOO_METIS_LIB}")
        list(APPEND MOO_EXTRA_LIBS ${MOO_METIS_LIB})
    else()
        message(FATAL_ERROR "METIS library not found")
    endif()
endif()

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

###### START HANDLE NON MSYS MAKEFILES ON WIN

if(WIN32)
    if(NOT CMAKE_GENERATOR STREQUAL "MSYS Makefiles")
        if(CMAKE_VERSION VERSION_LESS "3.21")
            message(FATAL_ERROR "Building MOO on Windows with Generator != 'MSYS Makefiles' requires CMake >= 3.21 for FindMsys support")
        endif()

        find_package(Msys) # setting 'REQUIRED', leads to failures somehow
        find_program(MSYS_BASH bash PATHS ${MSYS_INSTALL_PATH}/usr/bin REQUIRED)
        set(MOO_SHELL "${MSYS_BASH}")
        set(MOO_DOWNLOAD_SHELL "${MSYS_BASH}")
    else()
        set(MOO_SHELL "")
        set(MOO_DOWNLOAD_SHELL "sh")
    endif()
else()
    set(MOO_SHELL "sh")
    set(MOO_DOWNLOAD_SHELL "sh")
endif()

message(STATUS "Using MOO_SHELL: ${MSYS_BASH}")
message(STATUS "Using MOO_DOWNLOAD_SHELL: ${MOO_DOWNLOAD_SHELL}")

# Function: make_posix_path(VAR_NAME OUT_VAR)
#   VAR_NAME   : the name of the input variable (e.g., CMAKE_C_COMPILER)
#   OUT_VAR    : the name of the output variable (e.g., CMAKE_C_COMPILER_POSIX)
function(make_posix_path VAR_NAME OUT_VAR)
    set(_in_value "${${VAR_NAME}}")

    if (WIN32 AND (NOT CMAKE_GENERATOR STREQUAL "MSYS Makefiles"))
        set(CYGPATH "${MSYS_INSTALL_PATH}/usr/bin/cygpath.exe")

        # Convert to Posix-style path
        execute_process(
            COMMAND ${CYGPATH} -u "${_in_value}"
            OUTPUT_VARIABLE _posix_value
            OUTPUT_STRIP_TRAILING_WHITESPACE
        )

        # Export to the requested output variable
        set(${OUT_VAR} "${_posix_value}" PARENT_SCOPE)
    else()
        # On Linux/macOS or in MSYS just copy the original
        set(${OUT_VAR} "${_in_value}" PARENT_SCOPE)
    endif()
endfunction()

make_posix_path(CMAKE_C_COMPILER CMAKE_C_COMPILER_POSIX)
make_posix_path(CMAKE_CXX_COMPILER CMAKE_CXX_COMPILER_POSIX)
make_posix_path(CMAKE_Fortran_COMPILER CMAKE_Fortran_COMPILER_POSIX)
make_posix_path(MOO_LAPACK_LIB_STRING MOO_LAPACK_LIB_STRING_POSIX)
make_posix_path(CMAKE_CURRENT_SOURCE_DIR CMAKE_CURRENT_SOURCE_DIR_POSIX)
make_posix_path(CMAKE_CURRENT_BINARY_DIR CMAKE_CURRENT_BINARY_DIR_POSIX)

###### END HANDLE NON MSYS MAKEFILES ON WIN

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

if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    # we need to add this for gfortran / f95 as there are size mismatches between metis and mumps type sizes
    list(APPEND MOO_MUMPS_FORTRAN_FLAGS -fallow-argument-mismatch)
endif()

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

# TODO: add option to use binaries of mumps and ipopt that are already available
#       or fetched from MOO release - targets: Linux, MacOS, Windows (MSYS compatible)
# TODO: enable build w/o mumps or ipopt - alternative NLP solvers
# TODO: remove Ipopt examples / tutorial as they bloat the build

# for redistribution, we always call this manually as not all platforms can run this script
if(NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/third-party/mumps/MUMPS)
    message(STATUS "MUMPS sources not found. Trying to download MUMPS sources")
    message(STATUS "Downloading MUMPS sources requires patch")
    find_program(PATCH_EXECUTABLE patch)
    if (PATCH_EXECUTABLE)
        message(STATUS "patch found: ${PATCH_EXECUTABLE}")
    else()
        message(WARNING "patch not found! MUMPS download might fail.")
    endif()
    message(STATUS "Downloading mumps sources...")
    execute_process(COMMAND ${MOO_DOWNLOAD_SHELL} ./get.Mumps WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/third-party/mumps)
else()
    message(STATUS "MUMPS sources found")
endif()

ExternalProject_Add(moo_coinmumps
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/third-party/mumps
    BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/third-party/mumps
    CONFIGURE_COMMAND ${MOO_SHELL}
        ${CMAKE_CURRENT_SOURCE_DIR_POSIX}/third-party/mumps/configure
        CC=${CMAKE_C_COMPILER_POSIX}
        CXX=${CMAKE_CXX_COMPILER_POSIX}
        FC=${CMAKE_Fortran_COMPILER_POSIX}
        ADD_CFLAGS=-O3
        ADD_FCFLAGS=${MOO_MUMPS_FORTRAN_FLAGS_STRING}
        ADD_CXXFLAGS=-O3
        --with-lapack=${MOO_LAPACK_LIB_STRING_POSIX}
        --enable-shared=no
        --enable-static=yes
    BUILD_COMMAND ${MOO_SUBMAKE_CMD}
    INSTALL_COMMAND ""
    BUILD_ALWAYS 0
    BUILD_BYPRODUCTS ${CMAKE_CURRENT_BINARY_DIR}/third-party/mumps/.libs/libcoinmumps.a
)

add_library(moo::coinmumps STATIC IMPORTED GLOBAL)

set_target_properties(moo::coinmumps PROPERTIES
    IMPORTED_LOCATION ${CMAKE_CURRENT_BINARY_DIR}/third-party/mumps/.libs/libcoinmumps.a
    INTERFACE_INCLUDE_DIRECTORIES ${CMAKE_CURRENT_SOURCE_DIR}/third-party/mumps/MUMPS/include
)

add_dependencies(moo::coinmumps moo_coinmumps)

set(MUMPS_LFLAGS_POSIX "-L${CMAKE_CURRENT_BINARY_DIR_POSIX}/third-party/mumps/.libs")

set(MUMPS_CFLAGS_POSIX
    -DCOIN_USE_MUMPS_MPI_H
    -I${CMAKE_CURRENT_BINARY_DIR_POSIX}/third-party/mumps
    -I${CMAKE_CURRENT_SOURCE_DIR_POSIX}/third-party/mumps/MUMPS/include
    -I${CMAKE_CURRENT_SOURCE_DIR_POSIX}/third-party/mumps/MUMPS/libseq
)

string(JOIN " " MUMPS_CFLAGS_STRING_POSIX ${MUMPS_CFLAGS_POSIX})

ExternalProject_Add(moo_ipopt
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/third-party/ipopt
    BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/third-party/ipopt
    CONFIGURE_COMMAND ${MOO_SHELL}
    ${CMAKE_CURRENT_SOURCE_DIR_POSIX}/third-party/ipopt/configure
        CC=${CMAKE_C_COMPILE_POSIX}
        CXX=${CMAKE_CXX_COMPILER_POSIX}
        FC=${CMAKE_Fortran_COMPILER_POSIX}
        ADD_CFLAGS=-O3
        ADD_FCFLAGS=-O3
        ADD_CXXFLAGS=-O3
        --with-lapack=${MOO_LAPACK_LIB_STRING_POSIX}
        --without-pardiso
        --with-mumps
        --with-mumps-lflags=${MUMPS_LFLAGS_POSIX}
        --with-mumps-cflags=${MUMPS_CFLAGS_STRING_POSIX}
        --enable-shared=no
        --enable-static=yes
        --disable-sipopt
        --disable-java
    BUILD_COMMAND ${MOO_SUBMAKE_CMD}
    INSTALL_COMMAND ""
    BUILD_ALWAYS 0
    BUILD_BYPRODUCTS ${CMAKE_CURRENT_BINARY_DIR}/third-party/ipopt/src/.libs/libipopt.a
)

add_library(moo::ipopt STATIC IMPORTED GLOBAL)

set(IPOPT_INCLUDE_DIRS
    ${CMAKE_CURRENT_SOURCE_DIR}/third-party/ipopt/src/Interfaces
    ${CMAKE_CURRENT_SOURCE_DIR}/third-party/ipopt/src/Common
    ${CMAKE_CURRENT_SOURCE_DIR}/third-party/ipopt/src/Algorithm
    ${CMAKE_CURRENT_SOURCE_DIR}/third-party/ipopt/src/LinAlg
)

set_target_properties(moo::ipopt PROPERTIES
    IMPORTED_LOCATION ${CMAKE_CURRENT_BINARY_DIR}/third-party/ipopt/src/.libs/libipopt.a
    INTERFACE_INCLUDE_DIRECTORIES "${IPOPT_INCLUDE_DIRS}"
)

add_dependencies(moo::ipopt moo_ipopt)

###### END MUMPS + Ipopt

###### START RADAU5

if(MOO_WITH_RADAU)
    add_subdirectory(third-party/RADAU)
endif()

##### END RADAU5

target_include_directories(moo PUBLIC include)
target_include_directories(moo PUBLIC src)

target_include_directories(moo PRIVATE ${THIRD_PARTY_INSTALL_PREFIX}/include/coin-or)

target_link_directories(moo PRIVATE ${THIRD_PARTY_INSTALL_PREFIX}/lib)

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
target_link_libraries(moo PRIVATE moo::ipopt moo::coinmumps)
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

# Tests

if(MOO_TESTS)
    enable_testing()
    message(STATUS "CTest enabled")
endif()

if(MOO_WITH_GDOPT AND MOO_TESTS)
    add_executable(test_gdopt
        src/runner/entry.cpp src/runner/generated.c
    )

    target_include_directories(test_gdopt PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/include
        ${CMAKE_CURRENT_SOURCE_DIR}/src
    )
    target_link_libraries(test_gdopt PRIVATE moo)

    add_test(NAME test_gdopt COMMAND test_gdopt)
endif()
