cmake_minimum_required(VERSION 3.14)

if(NOT MOO_WITH_GDOPT)
    message(FATAL_ERROR "-DMOO_WITH_GDOPT=ON is required to build tests (-DMOO_TESTS=ON).")
endif()

if(NOT MOO_WITH_RADAU)
    message(FATAL_ERROR "-DMOO_WITH_RADAU=ON is required to build tests (-DMOO_TESTS=ON).")
endif()

add_executable(test_gdopt
    gdopt/entry.cpp
    gdopt/generated.c
)

target_include_directories(test_gdopt PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/gdopt
)

target_link_libraries(test_gdopt PRIVATE moo)

add_test(NAME test_gdopt COMMAND test_gdopt)
