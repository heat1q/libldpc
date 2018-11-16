#include "ldpc_types.h"
#include "scm_types.h"

void ldpc_code_t_st_setup(ldpc_code_t* code, const size_t ST_MAX_SIZE);
void lpdc_code_t_stopping_sets(ldpc_code_t* code);
size_t lpdc_code_t_erasure_decoding(ldpc_code_t* code, bits_t** in_bits, size_t** set);
void generate_submatrix(ldpc_code_t* code, ldpc_code_t* erasure_code, size_t* epsilon, const size_t epsilon_size);
