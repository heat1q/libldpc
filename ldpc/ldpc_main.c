#include "unistd.h"

#include "ldpc_types.h"
#include "functions.h"
#include "ldpc_cycles.h"
#include "ldpc_decoder.h"

int main(int argc, char* argv[])
{
    if(argc < 2)
    {
        printf("Code file needed!\n");
        exit(EXIT_FAILURE);
    }

    if (access(argv[1], F_OK) == -1)
    {
        printf("Cannot find code file %s!\n", argv[1]);
        exit(EXIT_FAILURE);
    }

    ldpc_code_t *code = calloc(1, sizeof(ldpc_code_t));

    read_ldpc_file(code, argv[1]);


    //girth_ldpc_code_t(code);
    //cycle_ldpc_code_t(code);

    printf("Test\n");
    double* llr_out = calloc(code->nc, sizeof(double));
    double* llr_in = calloc(code->nc, sizeof(double));
    for (size_t i = 0; i < code->nc; ++i)
        llr_in[i] = 0.5;

    ldpc_decode(*code, llr_in, &llr_out, 10000, 0);

    printf("[");
    for (size_t i = 0; i < code->nc; ++i)
    {
        printf("%.2f", llr_out[i]);
        if (i < code->nc - 1)
            printf(" ");
    }
    printf("]\n");

    destroy_ldpc_code_t(code);
	free(code);
    free(llr_out);
    free(llr_in);

    return 0;
}

