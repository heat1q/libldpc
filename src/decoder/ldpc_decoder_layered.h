#pragma once

#include "ldpc_decoder.h"

void layered_dec(ldpc_code_t* code, double* llr_in, double* l_c2v, double* l_c2v_sum, double* l_v2c, uint64_t* cn_subset, const uint64_t cn_size, double* f, double* b);
uint8_t layered_dec_setup(ldpc_decoder_lyr_t *dec, ldpc_code_t* code, char* clfile);
uint64_t ldpc_decode_layered(ldpc_decoder_lyr_t *dec, ldpc_code_t* code, double *llr_in, double* llr_out, const uint64_t MaxIter, const uint8_t early_termination);
void destroy_dec(ldpc_decoder_lyr_t *dec, ldpc_code_t *code);
