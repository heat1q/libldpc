#pragma once

#define MAX_FILENAME_LEN 256

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>

#ifdef INTELMKL
#include <mkl.h>
#endif


#ifndef M_PI
#define M_PI 3.14159265358979323846264338327
#endif

struct {
    double *pX;
    double *X;
    double *A;
    uint16_t M;
    uint16_t log2M;
} typedef cstll_t;

struct {
    uint64_t num_bits_in;
    uint64_t num_symbols_out;
    uint64_t* types;
    uint64_t num_types;
} typedef ccdm_t;

#ifdef SMDM
struct {
    double rate;
    size_t k;
    size_t n;
    double* pA;
    size_t k0;
    size_t n0;
} typedef dm_t;
#else
struct {
    double rate;
    size_t k;
    size_t n;
    double* pA;
    size_t* types;
    size_t num_types;
} typedef dm_t;
#endif


typedef uint8_t bits_t;
typedef uint16_t labels_t;
typedef uint32_t symbols_t;
