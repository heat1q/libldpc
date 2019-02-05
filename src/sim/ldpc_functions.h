#include <math.h>
#include <stdint.h>
#include "ldpc_types.h"
#include "scm_types.h"

#ifdef INTELMKL
void ldpc_encode(ldpc_sim_t sim, ldpc_code_t code, uint64_t* x, bits_t* c, VSLStreamStatePtr* rng_stream);
void ldpc_encode_all0(ldpc_sim_t sim, ldpc_code_t code, uint64_t* x, bits_t* c, VSLStreamStatePtr* rng_stream);
void ldpc_pas_encode(ldpc_sim_t sim, ldpc_code_t code, dm_t dm, uint64_t *x, bits_t* c, bits_t* data, uint8_t* symbols, VSLStreamStatePtr* rng_stream);
void ldpc_pas_encode_all0(ldpc_sim_t sim, ldpc_code_t code, dm_t dm, uint64_t *x, bits_t* c, bits_t* data, uint8_t* symbols, VSLStreamStatePtr* rng_stream);
#endif

void encode(ldpc_code_t code, bits_t* c);
void map_c_to_x(ldpc_sim_t sim, bits_t* c, size_t* x);
void ldpc_encode_all0_legacy(ldpc_sim_t sim, ldpc_code_t code, uint64_t* x, bits_t* c);

void calc_syndrome_c(ldpc_code_t code, bits_t* c, bits_t* s);
void calc_syndrome_llr(ldpc_code_t code, double* llr_out, bits_t* s);
size_t analyze_codeword(ldpc_code_t code, bits_t* c, size_t* one_pos);
bits_t is_codeword(ldpc_code_t code, bits_t* c);


void ldpc_pas_encode_legacy(ldpc_sim_t sim, ldpc_code_t code, dm_t dm, uint64_t *x, bits_t* c, bits_t* data, uint8_t* symbols);
void ldpc_pas_encode_all0_legacy(ldpc_sim_t sim, ldpc_code_t code, dm_t dm, uint64_t *x, bits_t* c, bits_t* data, uint8_t* symbols);

void setup_ldpc_sim(ldpc_sim_t *sim, ldpc_code_t* code, cstll_t* cstll, char *simfile, char *ldpcfile, char* mappingfile);
void setup_ldpc_sim_pas(ldpc_sim_t *sim, ldpc_code_t* code, cstll_t* cstll, dm_t *dm, char *simfile, char *ldpcfile, char* mappingfile);

int read_ldpc_file(ldpc_code_t* code, const char* filename);

void read_bit_mapping_file(char* filename, size_t** bit_mapping, size_t bits, size_t n);

void log_error(ldpc_sim_t sim, ldpc_code_t code, cstll_t cstll, bits_t *c, double *lout, size_t frame_num, double snr);

// void get_chat(ldpc_code_t code, double* llr_out, bits_t* chat);
// uint64_t cmp_codewords(bits_t* c1, bits_t* c2, size_t len);

void destroy_ldpc_code_t(ldpc_code_t* code);
void destroy_ldpc_sim_t(ldpc_sim_t *sim);

void print_ldpc_code_t(ldpc_code_t code);
void print_ldpc_sim_t(ldpc_sim_t sim);
