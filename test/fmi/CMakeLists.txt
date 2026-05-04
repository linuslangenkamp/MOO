add_executable(fmi_sanity_check sanity_check.c)

target_link_libraries(fmi_sanity_check PRIVATE fmi4c)

target_compile_definitions(fmi_sanity_check PRIVATE
    FMU_PATH="${CMAKE_CURRENT_SOURCE_DIR}/resources/BouncingBall.fmu"
)

add_test(NAME fmi_sanity_check COMMAND fmi_sanity_check)
