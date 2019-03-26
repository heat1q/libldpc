#pragma once

#include "ldpc_decoder.h"

uint8_t layered_dec_setup(ldpc_decoder_lyr_t *dec, ldpc_code_t* code, char* clfile);
uint64_t ldpc_decode_layered(ldpc_decoder_lyr_t *dec, ldpc_code_t* code, double *llr_in, double* llr_out, const uint64_t MaxIter, const uint8_t early_termination);
void destroy_dec(ldpc_decoder_lyr_t *dec, ldpc_code_t *code);
