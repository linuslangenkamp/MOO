cmake_minimum_required(VERSION 3.10)

option(WITH_RADAU5 "Build and link with 'RADAU5' library" TRUE)

if(WITH_RADAU5)
    message(STATUS "Library 'RADAU5' will be built from source and linked (WITH_RADAU5=TRUE)")

    project(radau5)

    enable_language(Fortran)

    find_package(LAPACK REQUIRED)

    add_library(radau5 STATIC
        radau5.f
        dc_lapack.f
    )

    target_compile_options(radau5 PUBLIC "-fPIC" "-O3")

    target_link_libraries(radau5 PUBLIC ${LAPACK_LIBRARIES} gfortran)

    target_include_directories(radau5 PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

else()
  message(STATUS "Library 'WITH_RADAU5' is disabled (WITH_RADAU5=FALSE)")
endif()

