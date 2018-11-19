#include "unistd.h"
#include "string.h"

#include "ldpc_types.h"
#include "functions.h"
#include "ldpc_cycles.h"
#include "ldpc_decoder.h"
#include "ldpc_stoppingsets.h"

int main(int argc, char* argv[])
{
    if(argc < 3)
    {
        printf("Supply code file and maximum size!\n");
        exit(EXIT_FAILURE);
    }

    if (access(argv[1], F_OK) == -1)
    {
        printf("Cannot find code file %s!\n", argv[1]);
        exit(EXIT_FAILURE);
    }

    ldpc_code_t *code = calloc(1, sizeof(ldpc_code_t));

    read_ldpc_file(code, argv[1]);

    lpdc_code_t_stopping_sets(code, (uint64_t) strtol(argv[2], NULL, 10));

    destroy_ldpc_code_t(code);
    free(code);

    return 0;
}

