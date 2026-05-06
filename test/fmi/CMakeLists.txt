# FMI sanity check
add_executable(fmi_sanity_check sanity_check.c)

target_link_libraries(fmi_sanity_check PRIVATE fmi4c)

target_compile_definitions(fmi_sanity_check PRIVATE
    FMU_PATH="${CMAKE_CURRENT_SOURCE_DIR}/resources/BouncingBall.fmu"
)

add_test(NAME fmi_sanity_check COMMAND fmi_sanity_check)

# fmi-ls-dae sanity check - references an FMU that is not redistributable,
# so the path is supplied via the SIMPLE_DAE_FMU environment variable.
if(DEFINED ENV{SIMPLE_DAE_FMU})
    add_executable(lsdae_sanity_check lsdae_sanity_check.c)
    target_link_libraries(lsdae_sanity_check PRIVATE fmi4c)
    target_compile_definitions(lsdae_sanity_check PRIVATE
        SIMPLE_DAE_FMU_PATH="$ENV{SIMPLE_DAE_FMU}"
    )
    add_test(NAME lsdae_sanity_check COMMAND lsdae_sanity_check)
else()
    message(WARNING "SIMPLE_DAE_FMU not set - skipping lsdae_sanity_check")
endif()

# FMI ls-dae test 1
add_executable(ls_dae_fourbar ls_dae_fourbar.cpp)

target_link_libraries(ls_dae_fourbar PRIVATE fmi4c moo)

target_compile_definitions(ls_dae_fourbar PRIVATE
    FOURBAR_FMU_PATH="$ENV{FOURBAR_FMU}"
)

add_test(NAME ls_dae_fourbar COMMAND ls_dae_fourbar)

# FMI ls-dae test 2
add_executable(ls_dae_simple_dae ls_dae_simple_dae.cpp)

target_link_libraries(ls_dae_simple_dae PRIVATE fmi4c moo)

target_compile_definitions(ls_dae_simple_dae PRIVATE
    SIMPLE_DAE_FMU_PATH="$ENV{SIMPLE_DAE_FMU}"
)

add_test(NAME ls_dae_simple_dae COMMAND ls_dae_simple_dae)
