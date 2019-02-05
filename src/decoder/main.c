#include <unistd.h>
#include <string.h>

#include "../function/ldpc_types.h"
#include "../function/ldpc_functions.h"
#include "ldpc_decoder_layered.h"

int main()
{
    ldpc_code_t *code = calloc(1, sizeof(ldpc_code_t));

    read_ldpc_file(code, "code6x6.txt");

    double *llrin = calloc(6, sizeof (double));
    double *llrout = calloc(6, sizeof (double));
    llrin[0] = 1;llrin[1] = 1;llrin[2] = 1;
    llrin[3] = 1;llrin[4] = 1;llrin[5] = 1;

    uint64_t vn_subset[6] = {0, 1, 2, 3, 4, 5};

    uint64_t cn_subset[3] = {0, 2, 4};
    ldpc_decode_layered(code, llrin, &llrout, 10, vn_subset, 6, cn_subset, 3);
    printf("Subset 1 :: \t");
    printVectorDouble(llrout, 6);

    cn_subset[0] = 1;cn_subset[0] = 3;cn_subset[0] = 5;
    ldpc_decode_layered(code, llrout, &llrout, 10, vn_subset, 6, cn_subset, 3);
    printf("Subset 2 :: \t");
    printVectorDouble(llrout, 6);

    ldpc_decode(*code, llrin, &llrout, 100, 0);
    printf("Full :: \t");
    printVectorDouble(llrout, 6);

    destroy_ldpc_code_t(code);
    free(code);
    free(llrin);
    free(llrout);

    return 0;
}
