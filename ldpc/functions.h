#include "ldpc_types.h"
#include "scm_types.h"

void read_ldpc_file(ldpc_code_t* code, char* filename);
void print_ldpc_code_t(ldpc_code_t code);
void destroy_ldpc_code_t(ldpc_code_t* code);

void calc_syndrome_c(ldpc_code_t code, bits_t* c, bits_t* s);
bits_t is_codeword(ldpc_code_t code, bits_t* c);

void printVector(void *input, const size_t k);
void printVectorDouble(double* x, const size_t k);
void printBits(bits_t* x, const size_t k);

void generic_codeword_search(ldpc_code_t* code, bits_t** bits, const size_t length, const size_t currentbit);
