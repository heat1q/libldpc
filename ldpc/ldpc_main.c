#include "unistd.h"

#include "ldpc_types.h"
#include "functions.h"
#include "ldpc_cycles.h"
#include "ldpc_decoder.h"
#include "ldpc_stoppingsets.h"

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

    //lpdc_code_t_stopping_sets(code)
    bits_t* in_bits = calloc(code->nc, sizeof(size_t));
    bits_t* out_bits = calloc(code->nc, sizeof(size_t));

    //generic_codeword_search(code, &in_bits, code->nc, 0);

    for (size_t i = 0; i < code->nc; ++i)
    {
        int bit;
        printf("Enter %lu. bit: ", i + 1);
        scanf("%i", &bit);
        in_bits[i] = (unsigned short)bit;
    }

    lpdc_code_t_erasure_decoding(code, in_bits, &out_bits);
    printBits(out_bits, code->nc);

    free(in_bits);
    free(out_bits);
    destroy_ldpc_code_t(code);
    free(code);

    return 0;
}

