#include <unistd.h>
#include <string.h>
#include <math.h>

#include "../function/ldpc_types.h"
#include "../function/ldpc_functions.h"
#include "ldpc_decoder_layered.h"

int main()
{
    ldpc_code_t *code = calloc(1, sizeof(ldpc_code_t));

    read_ldpc_file(code, "code9x9.txt");

    bits_t* bits = calloc(code->nc, sizeof (bits_t));
    double *llr_in = calloc(code->nc, sizeof (double));
    double *llr_out = calloc(code->nc, sizeof (double));
    for (size_t i = 0; i < code->nc; ++i)
    {
        llr_in[i] = (((double)(rand() % 100)) - 50.0)/2;
    }


    layered_dec_setup(code, "layer9x9.txt");
    ldpc_decode_layered(code, llr_in, llr_out, 5, 0);
    for (size_t i = 0; i < code->nc; ++i)
        bits[i] = (llr_out[i] <= 0);
    printf("Layered :: isCodeword=%i \t", is_codeword(*code, bits));
    printVectorDouble(llr_out, code->nc);


    ldpc_decode(*code, llr_in, llr_out, 5, 0);
    for (size_t i = 0; i < code->nc; ++i)
        bits[i] = (llr_out[i] <= 0);
    printf("Full :: isCodeword=%i \t", is_codeword(*code, bits));
    printVectorDouble(llr_out, code->nc);


    /*
    destroy_ldpc_code_t(code);
    free(code);
    free(llr_out);
    */

    return 0;
}
