cmake_minimum_required(VERSION 3.14)

add_executable(test_init
    init_objectives.cpp
)

target_link_libraries(test_init PRIVATE moo)

add_test(NAME test_init COMMAND test_init)

add_executable(test_init_hessian
    init_hessian.cpp
)

target_link_libraries(test_init_hessian PRIVATE moo)

add_test(NAME test_init_hessian COMMAND test_init_hessian)

add_executable(test_init_benchmark
    init_benchmark.cpp
)

target_link_libraries(test_init_benchmark PRIVATE moo)

add_test(NAME test_init_benchmark COMMAND test_init_benchmark)
