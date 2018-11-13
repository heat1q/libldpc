#include "ldpc_types.h"
#include "scm_types.h"

void lpdc_code_t_stopping_sets(ldpc_code_t* code);
void lpdc_code_t_erasure_decoding(ldpc_code_t* code, bits_t* in_bits, bits_t** out_bits);
void generate_submatrix(ldpc_code_t* code, ldpc_code_t& erasure_code, size_t* epsilon);
