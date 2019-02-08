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
    double *llr2 = calloc(6, sizeof (double));
    llr1[0] = -1;llr1[1] = 2;llr1[2] = -3;
    llr1[3] = 4;llr1[4] = -5;llr1[5] = 6;

    llr2[0] = -1;llr2[1] = 2;llr2[2] = -3;
    llr2[3] = 4;llr2[4] = -5;llr2[5] = 6;

    ldpc_decode_layered_init(code, &llr1);
    printf("Layerd :: \t");
    printVectorDouble(llr1, 6);

    ldpc_decode(*code, llr2, &llr2, 1, 0);
    printf("Full :: \t");
    printVectorDouble(llr2, 6);

    destroy_ldpc_code_t(code);
    free(code);
    free(llr1);
    free(llr2);

    return 0;
}
