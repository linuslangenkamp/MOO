add_executable(test_ad_graph_core graph_core.cpp)
target_include_directories(test_ad_graph_core PRIVATE ${CMAKE_SOURCE_DIR}/src/ad)
target_compile_definitions(test_ad_graph_core
    PRIVATE
        MOO_SOURCE_DIR="${CMAKE_SOURCE_DIR}"
        MOO_AD_SOURCE_DIR="${CMAKE_SOURCE_DIR}/src/ad")
target_link_libraries(test_ad_graph_core PRIVATE mooad)
add_test(NAME test_ad_graph_core COMMAND test_ad_graph_core)

add_executable(test_ad_function_core function_core.cpp)
target_include_directories(test_ad_function_core PRIVATE ${CMAKE_SOURCE_DIR}/src/ad)
target_compile_definitions(test_ad_function_core
    PRIVATE
        MOO_SOURCE_DIR="${CMAKE_SOURCE_DIR}"
        MOO_AD_SOURCE_DIR="${CMAKE_SOURCE_DIR}/src/ad")
target_link_libraries(test_ad_function_core PRIVATE mooad)
add_test(NAME test_ad_function_core COMMAND test_ad_function_core)

add_executable(test_ad_function_call function_call.cpp)
target_include_directories(test_ad_function_call PRIVATE ${CMAKE_SOURCE_DIR}/src/ad)
target_compile_definitions(test_ad_function_call
    PRIVATE
        MOO_SOURCE_DIR="${CMAKE_SOURCE_DIR}"
        MOO_AD_SOURCE_DIR="${CMAKE_SOURCE_DIR}/src/ad")
target_link_libraries(test_ad_function_call PRIVATE mooad)
add_test(NAME test_ad_function_call COMMAND test_ad_function_call)

add_executable(test_ad_function_transform function_transform.cpp)
target_include_directories(test_ad_function_transform PRIVATE ${CMAKE_SOURCE_DIR}/src/ad)
target_compile_definitions(test_ad_function_transform
    PRIVATE
        MOO_SOURCE_DIR="${CMAKE_SOURCE_DIR}"
        MOO_AD_SOURCE_DIR="${CMAKE_SOURCE_DIR}/src/ad")
target_link_libraries(test_ad_function_transform PRIVATE mooad)
add_test(NAME test_ad_function_transform COMMAND test_ad_function_transform)

add_executable(test_ad_forward forward.cpp)
target_include_directories(test_ad_forward PRIVATE ${CMAKE_SOURCE_DIR}/src/ad)
target_compile_definitions(test_ad_forward
    PRIVATE
        MOO_SOURCE_DIR="${CMAKE_SOURCE_DIR}"
        MOO_AD_SOURCE_DIR="${CMAKE_SOURCE_DIR}/src/ad")
target_link_libraries(test_ad_forward PRIVATE mooad)
add_test(NAME test_ad_forward COMMAND test_ad_forward)

add_executable(test_ad_reverse reverse.cpp)
target_include_directories(test_ad_reverse PRIVATE ${CMAKE_SOURCE_DIR}/src/ad)
target_compile_definitions(test_ad_reverse
    PRIVATE
        MOO_SOURCE_DIR="${CMAKE_SOURCE_DIR}"
        MOO_AD_SOURCE_DIR="${CMAKE_SOURCE_DIR}/src/ad")
target_link_libraries(test_ad_reverse PRIVATE mooad)
add_test(NAME test_ad_reverse COMMAND test_ad_reverse)

add_executable(test_ad_hvp_composition hvp_composition.cpp)
target_include_directories(test_ad_hvp_composition PRIVATE ${CMAKE_SOURCE_DIR}/src/ad)
target_compile_definitions(test_ad_hvp_composition
    PRIVATE
        MOO_SOURCE_DIR="${CMAKE_SOURCE_DIR}"
        MOO_AD_SOURCE_DIR="${CMAKE_SOURCE_DIR}/src/ad")
target_link_libraries(test_ad_hvp_composition PRIVATE mooad)
add_test(NAME test_ad_hvp_composition COMMAND test_ad_hvp_composition)

add_executable(test_ad_sparsity sparsity.cpp)
target_include_directories(test_ad_sparsity PRIVATE ${CMAKE_SOURCE_DIR}/src/ad)
target_compile_definitions(test_ad_sparsity
    PRIVATE
        MOO_SOURCE_DIR="${CMAKE_SOURCE_DIR}"
        MOO_AD_SOURCE_DIR="${CMAKE_SOURCE_DIR}/src/ad")
target_link_libraries(test_ad_sparsity PRIVATE mooad)
add_test(NAME test_ad_sparsity COMMAND test_ad_sparsity)

add_executable(test_ad_mapping_core mapping_core.cpp)
target_include_directories(test_ad_mapping_core PRIVATE ${CMAKE_SOURCE_DIR}/src/ad)
target_compile_definitions(test_ad_mapping_core
    PRIVATE
        MOO_SOURCE_DIR="${CMAKE_SOURCE_DIR}"
        MOO_AD_SOURCE_DIR="${CMAKE_SOURCE_DIR}/src/ad")
target_link_libraries(test_ad_mapping_core PRIVATE mooad)
add_test(NAME test_ad_mapping_core COMMAND test_ad_mapping_core)

add_executable(test_ad_mapped_derivatives mapped_derivatives.cpp)
target_include_directories(test_ad_mapped_derivatives PRIVATE ${CMAKE_SOURCE_DIR}/src/ad)
target_compile_definitions(test_ad_mapped_derivatives
    PRIVATE
        MOO_SOURCE_DIR="${CMAKE_SOURCE_DIR}"
        MOO_AD_SOURCE_DIR="${CMAKE_SOURCE_DIR}/src/ad")
target_link_libraries(test_ad_mapped_derivatives PRIVATE mooad)
add_test(NAME test_ad_mapped_derivatives COMMAND test_ad_mapped_derivatives)

add_executable(test_ad_mapped_sparsity mapped_sparsity.cpp)
target_include_directories(test_ad_mapped_sparsity PRIVATE ${CMAKE_SOURCE_DIR}/src/ad)
target_compile_definitions(test_ad_mapped_sparsity
    PRIVATE
        MOO_SOURCE_DIR="${CMAKE_SOURCE_DIR}"
        MOO_AD_SOURCE_DIR="${CMAKE_SOURCE_DIR}/src/ad")
target_link_libraries(test_ad_mapped_sparsity PRIVATE mooad)
add_test(NAME test_ad_mapped_sparsity COMMAND test_ad_mapped_sparsity)

add_executable(test_ad_global_assembly global_assembly.cpp)
target_include_directories(test_ad_global_assembly PRIVATE ${CMAKE_SOURCE_DIR}/src/ad)
target_compile_definitions(test_ad_global_assembly
    PRIVATE
        MOO_SOURCE_DIR="${CMAKE_SOURCE_DIR}"
        MOO_AD_SOURCE_DIR="${CMAKE_SOURCE_DIR}/src/ad")
target_link_libraries(test_ad_global_assembly PRIVATE mooad)
add_test(NAME test_ad_global_assembly COMMAND test_ad_global_assembly)

add_executable(test_ad_eval eval.cpp)
target_include_directories(test_ad_eval PRIVATE ${CMAKE_SOURCE_DIR}/src/ad)
target_compile_definitions(test_ad_eval
    PRIVATE
        MOO_SOURCE_DIR="${CMAKE_SOURCE_DIR}"
        MOO_AD_SOURCE_DIR="${CMAKE_SOURCE_DIR}/src/ad")
target_link_libraries(test_ad_eval PRIVATE mooad)
add_test(NAME test_ad_eval COMMAND test_ad_eval)

add_executable(test_ad_codegen codegen.cpp)
target_include_directories(test_ad_codegen PRIVATE ${CMAKE_SOURCE_DIR}/src/ad)
target_compile_definitions(test_ad_codegen
    PRIVATE
        MOO_SOURCE_DIR="${CMAKE_SOURCE_DIR}"
        MOO_AD_SOURCE_DIR="${CMAKE_SOURCE_DIR}/src/ad"
        MOO_BINARY_DIR="${CMAKE_BINARY_DIR}"
        MOO_C_COMPILER="${CMAKE_C_COMPILER}")
target_link_libraries(test_ad_codegen PRIVATE mooad)
add_test(NAME test_ad_codegen COMMAND test_ad_codegen)

add_executable(test_ad_simplification simplification.cpp)
target_include_directories(test_ad_simplification PRIVATE ${CMAKE_SOURCE_DIR}/src/ad)
target_compile_definitions(test_ad_simplification
    PRIVATE
        MOO_SOURCE_DIR="${CMAKE_SOURCE_DIR}"
        MOO_AD_SOURCE_DIR="${CMAKE_SOURCE_DIR}/src/ad")
target_link_libraries(test_ad_simplification PRIVATE mooad)
add_test(NAME test_ad_simplification COMMAND test_ad_simplification)

add_executable(test_ad_elementary_ops elementary_ops.cpp)
target_include_directories(test_ad_elementary_ops PRIVATE ${CMAKE_SOURCE_DIR}/src/ad)
target_compile_definitions(test_ad_elementary_ops
    PRIVATE
        MOO_SOURCE_DIR="${CMAKE_SOURCE_DIR}"
        MOO_AD_SOURCE_DIR="${CMAKE_SOURCE_DIR}/src/ad"
        MOO_BINARY_DIR="${CMAKE_BINARY_DIR}"
        MOO_C_COMPILER="${CMAKE_C_COMPILER}")
target_link_libraries(test_ad_elementary_ops PRIVATE mooad)
add_test(NAME test_ad_elementary_ops COMMAND test_ad_elementary_ops)

add_executable(test_ad_multi_output_function multi_output_function.cpp)
target_include_directories(test_ad_multi_output_function PRIVATE ${CMAKE_SOURCE_DIR}/src/ad)
target_compile_definitions(test_ad_multi_output_function
    PRIVATE
        MOO_SOURCE_DIR="${CMAKE_SOURCE_DIR}"
        MOO_AD_SOURCE_DIR="${CMAKE_SOURCE_DIR}/src/ad"
        MOO_BINARY_DIR="${CMAKE_BINARY_DIR}"
        MOO_C_COMPILER="${CMAKE_C_COMPILER}")
target_link_libraries(test_ad_multi_output_function PRIVATE mooad)
add_test(NAME test_ad_multi_output_function COMMAND test_ad_multi_output_function)

add_executable(test_ad_matrix_views matrix_views.cpp)
target_include_directories(test_ad_matrix_views PRIVATE ${CMAKE_SOURCE_DIR}/src/ad)
target_compile_definitions(test_ad_matrix_views
    PRIVATE
        MOO_SOURCE_DIR="${CMAKE_SOURCE_DIR}"
        MOO_AD_SOURCE_DIR="${CMAKE_SOURCE_DIR}/src/ad"
        MOO_BINARY_DIR="${CMAKE_BINARY_DIR}"
        MOO_C_COMPILER="${CMAKE_C_COMPILER}")
target_link_libraries(test_ad_matrix_views PRIVATE mooad)
add_test(NAME test_ad_matrix_views COMMAND test_ad_matrix_views)

add_executable(test_ad_symbolic_matrix_ops symbolic_matrix_ops.cpp)
target_include_directories(test_ad_symbolic_matrix_ops PRIVATE ${CMAKE_SOURCE_DIR}/src/ad)
target_compile_definitions(test_ad_symbolic_matrix_ops
    PRIVATE
        MOO_SOURCE_DIR="${CMAKE_SOURCE_DIR}"
        MOO_AD_SOURCE_DIR="${CMAKE_SOURCE_DIR}/src/ad"
        MOO_BINARY_DIR="${CMAKE_BINARY_DIR}"
        MOO_C_COMPILER="${CMAKE_C_COMPILER}")
target_link_libraries(test_ad_symbolic_matrix_ops PRIVATE mooad)
add_test(NAME test_ad_symbolic_matrix_ops COMMAND test_ad_symbolic_matrix_ops)

add_executable(test_ad_mapaccum map_accum.cpp)
target_include_directories(test_ad_mapaccum PRIVATE ${CMAKE_SOURCE_DIR}/src/ad)
target_compile_definitions(test_ad_mapaccum
    PRIVATE
        MOO_SOURCE_DIR="${CMAKE_SOURCE_DIR}"
        MOO_AD_SOURCE_DIR="${CMAKE_SOURCE_DIR}/src/ad"
        MOO_BINARY_DIR="${CMAKE_BINARY_DIR}"
        MOO_C_COMPILER="${CMAKE_C_COMPILER}")
target_link_libraries(test_ad_mapaccum PRIVATE mooad)
add_test(NAME test_ad_mapaccum COMMAND test_ad_mapaccum)

add_executable(test_ad_linear_solve linear_solve.cpp)
target_include_directories(test_ad_linear_solve PRIVATE ${CMAKE_SOURCE_DIR}/src/ad)
target_compile_definitions(test_ad_linear_solve
    PRIVATE
        MOO_SOURCE_DIR="${CMAKE_SOURCE_DIR}"
        MOO_AD_SOURCE_DIR="${CMAKE_SOURCE_DIR}/src/ad"
        MOO_BINARY_DIR="${CMAKE_BINARY_DIR}"
        MOO_C_COMPILER="${CMAKE_C_COMPILER}"
        MOO_CXX_COMPILER="${CMAKE_CXX_COMPILER}"
        MOO_AD_LIBRARY="$<TARGET_FILE:mooad>"
        MOO_AD_LIBRARY_DIR="$<TARGET_FILE_DIR:mooad>")
target_link_libraries(test_ad_linear_solve PRIVATE mooad)
add_test(NAME test_ad_linear_solve COMMAND test_ad_linear_solve)
