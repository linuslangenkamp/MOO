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

# FMI ls-dae test
add_executable(fmi_ls_dae ls_dae.cpp)

target_link_libraries(fmi_ls_dae PRIVATE fmi4c moo)

target_compile_definitions(fmi_ls_dae PRIVATE
    SIMPLE_DAE_FMU_PATH="$ENV{SIMPLE_DAE_FMU}"
)

add_test(NAME fmi_ls_dae COMMAND fmi_ls_dae)
