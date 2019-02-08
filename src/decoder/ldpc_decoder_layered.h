#include "ldpc_decoder.h"

void ldpc_decode_layered(ldpc_code_t* code, double** llr, uint64_t *cn_subset, const uint64_t cn_size, const size_t SubIter);
uint64_t ldpc_decode_layered_init(ldpc_code_t* code, double** llr, const uint64_t MaxIter);
