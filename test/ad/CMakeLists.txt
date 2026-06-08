add_executable(test_ad_sparsity sparsity.cpp)
target_include_directories(test_ad_sparsity PRIVATE ${CMAKE_SOURCE_DIR}/src/ad)
target_link_libraries(test_ad_sparsity PRIVATE mooad)
add_test(NAME test_ad_sparsity COMMAND test_ad_sparsity)

add_executable(test_ad_sparse_derivatives sparse_derivatives.cpp)
target_include_directories(test_ad_sparse_derivatives PRIVATE ${CMAKE_SOURCE_DIR}/src/ad)
target_link_libraries(test_ad_sparse_derivatives PRIVATE mooad)
add_test(NAME test_ad_sparse_derivatives COMMAND test_ad_sparse_derivatives)

add_executable(test_ad_vector_expressions vector_expressions.cpp)
target_include_directories(test_ad_vector_expressions PRIVATE ${CMAKE_SOURCE_DIR}/src/ad)
target_link_libraries(test_ad_vector_expressions PRIVATE mooad)
add_test(NAME test_ad_vector_expressions COMMAND test_ad_vector_expressions)

if(MOO_WITH_PYTHON_BINDINGS)
    add_test(NAME test_ad_python_bindings
             COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/python_bindings.py)
    set_tests_properties(test_ad_python_bindings PROPERTIES
        ENVIRONMENT "PYTHONPATH=${CMAKE_SOURCE_DIR}/python:${CMAKE_BINARY_DIR}/python")
endif()
