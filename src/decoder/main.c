#include <unistd.h>
#include <string.h>

#include "../function/ldpc_types.h"
#include "../function/ldpc_functions.h"
#include "ldpc_decoder_layered.h"

int main()
{
    ldpc_code_t *code = calloc(1, sizeof(ldpc_code_t));

    read_ldpc_file(code, "code6x6.txt");

    double *llr1 = calloc(6, sizeof (double));
    double *llr_out = calloc(6, sizeof (double));
    llr1[0] = -1;llr1[1] = 2;llr1[2] = -3;
    llr1[3] = 4;llr1[4] = -5;llr1[5] = 6;

    ldpc_decode_layered_init(code, llr1, 100);

    for(int i=0; i<code->nnz; ++i)
        llr_out[i] = 0.0;
    ldpc_decode(*code, llr1, &llr_out, 100, 0);
    printf("Full :: \t");
    printVectorDouble(llr_out, 6);

    /*
    destroy_ldpc_code_t(code);
    free(code);
    free(llr1);
    free(llr_out);
    */

    return 0;
}
