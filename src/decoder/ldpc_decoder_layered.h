#include "ldpc_decoder.h"

uint64_t ldpc_decode_layered(ldpc_code_t* code, double* llr_in, double** llr_out, uint64_t max_iter, uint64_t *vn_subset, const uint64_t vn_size, uint64_t *cn_subset, const uint64_t cn_size);
