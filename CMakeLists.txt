cmake_minimum_required(VERSION 3.14)
project(moo LANGUAGES C CXX Fortran)

set(CMAKE_CXX_STANDARD 17)

### Build command
# cmake .. -DMOO_BUILD_TYPE=Release -DCMAKE_BUILD_TYPE=Release -DMOO_WITH_GDOPT_EXAMPLE=1 -DMUMPS_USE_METIS=ON / OFF

### Options

option(MOO_WITH_RADAU "Enable RADAU build and linking" ON)
option(MOO_WITH_GDOPT "Build MOO with GDOPT support" ON)
option(MOO_WITH_GDOPT_EXAMPLE "Build MOO with GDOPT example executable" OFF)

##### START setup

message(STATUS "C Compiler: ${CMAKE_C_COMPILER}")
message(STATUS "C++ Compiler: ${CMAKE_CXX_COMPILER}")
message(STATUS "Fortran Compiler: ${CMAKE_Fortran_COMPILER}")

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    set(IS_GCC TRUE)
else()
    set(IS_GCC FALSE)
endif()

##### END setup

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

# dmumps, mumps_common, seq static libs
target_link_directories(moo PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/bin)

###### END MUMPS + Ipopt

###### START RADAU5

if(MOO_WITH_RADAU)
    add_subdirectory(third-party/RADAU)
endif()

##### END RADAU5

target_include_directories(moo PUBLIC include)
target_include_directories(moo PUBLIC src)

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

if(MOO_WITH_RADAU)
    target_link_libraries(moo PRIVATE radau)
endif()

# hide symbols
if(WIN32)
    target_compile_definitions(moo PRIVATE MOO_DLL_EXPORT)
elseif (${CMAKE_SYSTEM_NAME} STREQUAL "Darwin")
    target_link_libraries(moo PRIVATE dl)
    target_link_directories(moo PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/bin)
    target_link_options(moo PRIVATE -Wl,-hidden-ldmumps -Wl,-hidden-lipopt)
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
