#include <unistd.h>
#include <string.h>
#include <math.h>

#include "../function/ldpc_types.h"
#include "../function/ldpc_functions.h"
#include "ldpc_decoder_layered.h"

int main()
{
    ldpc_code_t code;

    read_ldpc_file(&code, "../code/test_code/code_rand_proto_3x6_400_4.txt");

    bits_t* bits = calloc(code.nc, sizeof (bits_t));
    double *llr_in = calloc(code.nc, sizeof (double));
    double *llr_out = calloc(code.nc, sizeof (double));
    for (size_t i = 0; i < code.nc; ++i)
    {
        llr_in[i] = (((double)(rand() % 100)) - 50.0)/2;
    }


    layered_dec_setup(&code, "../code/test_code/layer_rand_proto_3x6_400_4.txt");
    ldpc_decode_layered(&code, llr_in, llr_out, 200, 0);
    for (size_t i = 0; i < code.nc; ++i)
        bits[i] = (llr_out[i] <= 0);
    printf("Layered :: ");
    printVectorDouble(llr_out, code.nc);
    printf("isCodeword=%i :: ", is_codeword(code, bits));
    printBits(bits, code.nc);


    ldpc_decode(code, llr_in, llr_out, 200, 0);
    for (size_t i = 0; i < code.nc; ++i)
        bits[i] = (llr_out[i] <= 0);
    printf("Full :: ");
    printVectorDouble(llr_out, code.nc);
    printf("isCodeword=%i :: ", is_codeword(code, bits));
    printBits(bits, code.nc);


    /*
    destroy_ldpc_code_t(code);
    free(code);
    free(llr_out);
    */

    return 0;
}
