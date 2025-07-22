file(GLOB SRCS *.c *.h)

add_executable(include_test_cmk ${SRCS})

target_include_directories(include_test_cmk PRIVATE
  /home/linus/Projects/OpenModelica/build_cmake/install_cmake/include/omc/c
  /home/linus/Projects/OpenModelica/build_cmake/install_cmake/include/omc
  ${CMAKE_CURRENT_SOURCE_DIR}
)

target_compile_definitions(include_test_cmk PRIVATE
  OM_HAVE_PTHREADS
  OMC_DLL_MAIN_DEFINE
  OPENMODELICA_XML_FROM_FILE_AT_RUNTIME
  OMC_MODEL_PREFIX=include_test
  OMC_NUM_MIXED_SYSTEMS=0
  OMC_NUM_LINEAR_SYSTEMS=0
  OMC_NUM_NONLINEAR_SYSTEMS=0
  OMC_NDELAY_EXPRESSIONS=0
  OMC_NVAR_STRING=0
)

target_compile_options(include_test_cmk PRIVATE -O2 -fPIC)

target_link_directories(include_test_cmk PRIVATE
  /home/linus/Projects/OpenModelica/build_cmake/install_cmake/lib/x86_64-linux-gnu/omc
)

target_link_libraries(include_test_cmk PRIVATE
  MOO
  SimulationRuntimeC
  OpenModelicaRuntimeC
  omcgc
  zlib
  lapack
  blas
  m
  dl
  pthread
  gfortran
  stdc++
)

set_target_properties(include_test_cmk PROPERTIES
  INSTALL_RPATH "/home/linus/Projects/OpenModelica/build_cmake/install_cmake/lib"
)