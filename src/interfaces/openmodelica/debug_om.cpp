#include "debug_om.h"

void print_real_var_names(DATA* data) {
    for (long idx = 0; idx < data->modelData->nVariablesReal; idx++)
        printf("%s\n", data->modelData->realVarsData[idx].info.name);
}

void print_real_var_names_values(DATA* data) {
    printf("Time: %f\n", data->localData[0]->timeValue);
    for (long idx = 0; idx < data->modelData->nVariablesReal; idx++)
        printf("%s: %f\n", data->modelData->realVarsData[idx].info.name, data->localData[0]->realVars[idx]);
}

void print_jacobian_sparsity(const JACOBIAN* jac, bool print_pattern, const char* name = nullptr) {
    if (!jac || !jac->sparsePattern) {
        printf("Invalid JACOBIAN or missing SPARSE_PATTERN\n");
        return;
    }

    const SPARSE_PATTERN* sp = jac->sparsePattern;
    const unsigned int nRows = (unsigned int)jac->sizeRows;
    const unsigned int nCols = (unsigned int)jac->sizeCols;

    printf("\n=== JACOBIAN SPARSITY INFO ===\n");
    if (name) {
        printf("Name: %s\n", name);
    }
    printf("Jacobian: %u rows x %u cols\n", nRows, nCols);
    printf("numberOfNonZeros: %u\n", sp->numberOfNonZeros);
    printf("sizeofIndex:      %u\n", sp->sizeofIndex);
    printf("maxColors:        %u\n", sp->maxColors);

    printf("leadindex: ");
    for (unsigned int i = 0; i <= nCols; ++i) {
        printf("%u ", sp->leadindex ? sp->leadindex[i] : 0);
    }
    printf("\n");

    printf("index:     ");
    for (unsigned int i = 0; i < sp->sizeofIndex; ++i) {
        printf("%u ", sp->index ? sp->index[i] : 0);
    }
    printf("\n");

    printf("colorCols: ");
    for (unsigned int i = 0; i < sp->maxColors; ++i) {
        printf("%u ", sp->colorCols ? sp->colorCols[i] : 0);
    }
    printf("\n");

    if (!print_pattern) {
        printf("===============================\n");
        return;
    }

    printf("\n=== JACOBIAN SPARSITY PLOT ===\n");
printf("      "); for (unsigned int col = 0; col < nCols; ++col) printf("%u", col % 10); printf("\n");
    for (unsigned int row = 0; row < nRows; ++row) {
        printf("%4u: ", row); // row index
        for (unsigned int col = 0; col < nCols; ++col) {
            unsigned int col_start = sp->leadindex[col];
            unsigned int col_end = sp->leadindex[col + 1];
            bool found = false;
            for (unsigned int i = col_start; i < col_end; ++i) {
                if (sp->index[i] == row) {
                    found = true;
                    break;
                }
            }
            printf("%c", found ? '*' : ' ');
        }
        printf("\n");
    }

    printf("================================\n");
}
