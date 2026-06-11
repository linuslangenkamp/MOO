function(add_ad_test name)
    set(extra_defs ${ARGN})
    add_executable(test_ad_${name} ${name}.cpp)
    target_include_directories(test_ad_${name} PRIVATE ${CMAKE_SOURCE_DIR}/src/ad)
    target_compile_definitions(test_ad_${name}
        PRIVATE
            MOO_SOURCE_DIR="${CMAKE_SOURCE_DIR}"
            MOO_AD_SOURCE_DIR="${CMAKE_SOURCE_DIR}/src/ad"
            ${extra_defs})
    target_link_libraries(test_ad_${name} PRIVATE mooad)
    add_test(NAME test_ad_${name} COMMAND test_ad_${name})
endfunction()

set(CODEGEN_DEFS
    MOO_BINARY_DIR="${CMAKE_BINARY_DIR}"
    MOO_C_COMPILER="${CMAKE_C_COMPILER}")

set(LINEAR_SOLVE_DEFS
    ${CODEGEN_DEFS}
    MOO_CXX_COMPILER="${CMAKE_CXX_COMPILER}"
    MOO_AD_LIBRARY="$<TARGET_FILE:mooad>"
    MOO_AD_LIBRARY_DIR="$<TARGET_FILE_DIR:mooad>")

# basic functionality tests
add_ad_test(graph_core)
add_ad_test(function_core)
add_ad_test(function_call)
add_ad_test(function_transform)
add_ad_test(forward)
add_ad_test(reverse)
add_ad_test(hvp_composition)
add_ad_test(sparsity)
add_ad_test(mapping_core)
add_ad_test(mapped_derivatives)
add_ad_test(mapped_sparsity)
add_ad_test(global_assembly)
add_ad_test(eval)
add_ad_test(simplification)

# tests requiring codegen definitions
add_ad_test(codegen              ${CODEGEN_DEFS})
add_ad_test(elementary_ops       ${CODEGEN_DEFS})
add_ad_test(multi_output_function ${CODEGEN_DEFS})
add_ad_test(matrix_views         ${CODEGEN_DEFS})
add_ad_test(symbolic_matrix_ops  ${CODEGEN_DEFS})
add_ad_test(symbolic_sparse_matrix ${CODEGEN_DEFS})
add_ad_test(map_accum            ${CODEGEN_DEFS})

# linear solve (needs extra defs on top of codegen)
add_ad_test(linear_solve         ${LINEAR_SOLVE_DEFS})
