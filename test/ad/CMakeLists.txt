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
