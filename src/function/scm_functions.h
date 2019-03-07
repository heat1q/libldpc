#pragma once

#include "scm_types.h"

#define MAX_LLR 9999.9
#define MIN_LLR -9999.9

void calc_llrs(double y, cstll_t cstll, double sigma2, uint16_t* labels, double* llrs_out);
void calc_llrs_ml(double y, bits_t* bits, uint8_t level, cstll_t cstll, double sigma2, uint16_t* labels, double* llr_out);

uint8_t check_type(symbols_t* symbols, size_t len, uint64_t* type, size_t type_len);
uint64_t sum(uint64_t *summands, size_t len);
double fsum(double *summands, size_t len);

#ifdef INTELMKL
void simulate_awgn(cstll_t cstll, uint64_t *x, double *y, double sigma2, uint64_t nc, VSLStreamStatePtr* rng_stream);
#endif

double randn();
double simulate_awgn_legacy(cstll_t cstll, uint64_t *x, double *y, double sigma2, uint64_t nc);
void simulate_bec(uint64_t *x, double* y, double epsilon, uint64_t nc);

int8_t sign(double a);

void reverse_labels(uint16_t* labels, uint16_t* labels_rev, uint16_t M);

uint64_t ipow(uint64_t base, uint64_t expo);

void shuffle(size_t* idx, size_t len);
void create_nb_interleaver(size_t* pi, size_t permute_len, size_t len);
void interleaver(bits_t* input, bits_t* output, size_t* pi, size_t len);

void init_double_array(double *arr, size_t len, double val);

void destroy_cstll_t(cstll_t* cstll);
void destroy_dm_t(dm_t* dm);

void dec2bin(uint64_t val, uint8_t pad);
void print_cstll_t(cstll_t cstll);
void print_dm_t(dm_t dm);
