#include "../function/ldpc_types.h"

void cycle_ldpc_code_t(ldpc_code_t* code);
void girth_ldpc_code_t(ldpc_code_t* code);

void ex_msg_t_add(uint64_t* result, uint64_t* x, const size_t length);
void ex_msg_t_sum(uint64_t* result, uint64_t* data, const size_t length, const size_t j);
