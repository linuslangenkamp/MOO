cmake_minimum_required(VERSION 3.14)

add_executable(test_gdop
    entry.cpp
    generated.c
    sanity_check.c
    full_functionality_check.c
)

target_include_directories(test_gdop PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(test_gdop PRIVATE moo)

add_test(NAME test_gdop COMMAND test_gdop)
set_tests_properties(test_gdop PROPERTIES WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
