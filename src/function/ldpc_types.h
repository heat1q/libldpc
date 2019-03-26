#pragma once

#define _POSIX_C_SOURCE 199309L

#include "scm_types.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>
#include <time.h>
#include <math.h>

struct {
    uint64_t nc; /* blocklength of code */
    uint64_t kc; /* information length of code */
    uint64_t mc; /* number of parity check equations */
    uint64_t nct; /* number of transmitted code bits */
    uint64_t kct; /* number of transmitted information bits */
    uint64_t mct; /* number of transmitted parity check bits */
    uint64_t nnz; /* # of non-zero entries in H, i.e., # of edges in the Tanner graph */
    size_t* puncture; /* array pf punctured bit indices */
    size_t num_puncture; /* number of punctured bits */
    size_t num_puncture_sys; /* number of punctured bits in systematic part */
    size_t num_puncture_par; /* number of punctured bits in parity part */
    size_t* shorten; /* array of shortened bit indices */
    size_t num_shorten; /* number of shortened bits */
    uint64_t max_dc; /* maximum check node degree */
    uint8_t** genmat;
    size_t M;
    size_t N;
    size_t L;
    size_t mu;
    size_t window;
    size_t* cw; /* denotes the check weight of each check node, i.e., # of connected VN; dimensions cw[mc] */
    size_t* vw; /* denotes the variable weight, i.e., # of connected CN; dimensions vw[nc] */
    size_t** cn; /* denotes the check neighbors, i.e. connected VN, for each check node as index in c/r; dimensions cn[mc][cw[i]] */
    size_t** vn; /* denotes the var neighbors, i.e., connected CN, for each variable node as index in c/r; dimensions vn[nc][vw[i]] */
    size_t* r; /* non zero row indices; length nnz */
    size_t* c; /* non zero check indices; length nnz */
    uint64_t girth;
    size_t st_max_size;
    size_t* stw;
    char** st;
    //size_t* dw;
    //size_t** ds;
    uint64_t nl; //number of layers
    uint64_t* lw; //layer weight
    uint64_t** layers;

} typedef ldpc_code_t;

struct {
    uint64_t n;
    uint16_t M;
    uint16_t bits;
    uint64_t max_frames;
    uint64_t min_fec;
    uint64_t bp_iter;
    double* snrs;
    uint16_t* labels;
    uint16_t* labels_rev;
    uint8_t decoder_terminate_early;
    double SE;
    size_t num_snrs;
    char logfile[MAX_FILENAME_LEN];
    size_t** bit_mapper;
    size_t* bits_pos;
} typedef ldpc_sim_t;

struct {
    double* lsum;
    double* lsum_tmp;
    bits_t* c_out;

    double** l_c2v;
    double** l_v2c;
    double** f;
    double** b;
} typedef ldpc_decoder_lyr_t;
