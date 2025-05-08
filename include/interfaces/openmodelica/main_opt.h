#ifndef OPT_MAIN_H
#define OPT_MAIN_H

#include "simulation_data.h"

#ifdef __cplusplus
extern "C" {
#endif

int _main_OptimitationRuntime(int argc, char**argv, DATA *data, threadData_t *threadData);

#ifdef __cplusplus
};
#endif

#endif  // OPT_MAIN_H
