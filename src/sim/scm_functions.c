#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#ifdef INTELMKL
#include <mkl_vsl.h>
#include <mkl.h>
#endif
#include "scm_functions.h"


void calc_llrs(double y, cstll_t cstll, double sigma2, uint16_t* labels, double* llrs_out) {

    double tmp0, tmp1;

    for(size_t i = 0; i < cstll.log2M; i++) {
        tmp0 = 0.0;
        tmp1 = 0.0;
        for(size_t j = 0; j < cstll.M; j++) {
            if(labels[j] & (1 << (cstll.log2M-1-i))) {
                tmp1 += exp(-(y-cstll.X[j])*(y-cstll.X[j])/(2*sigma2)) * cstll.pX[j];
            } else {
                tmp0 += exp(-(y-cstll.X[j])*(y-cstll.X[j])/(2*sigma2)) * cstll.pX[j];
            }
        }
        double val = log(tmp0/tmp1);
        // check usually required when PAS is used with large constellations
        // and severely shaped distributions
        if(isinf(val) == +1) {
            llrs_out[i] = MAX_LLR;
        } else if(isinf(val) == -1) {
            llrs_out[i] = MIN_LLR;
        } else {
            llrs_out[i] = val;
        }
    }
}

void calc_llrs_ml(double y, bits_t* bits, uint8_t level, cstll_t cstll, double sigma2, uint16_t* labels, double* llr_out) {
    double tmp0, tmp1;
    int16_t mask0, mask1;

    tmp0 = tmp1 = 0;
    mask0 = mask1 = 0;

    /* create bit mask */
    for(size_t i = 0; i < level; i++) {
        mask0 += bits[i] << (level-i);
    }
    mask1 = mask0;
    mask1 += 1;

    for(size_t j = 0; j < cstll.M; j++) {
        if((labels[j] >> (cstll.log2M-level-1)) == mask1) {
            tmp1 += exp(-(y-cstll.X[j])*(y-cstll.X[j])/(2*sigma2)) * cstll.pX[j];
        } else if((labels[j] >> (cstll.log2M-level-1)) == mask0) {
            tmp0 += exp(-(y-cstll.X[j])*(y-cstll.X[j])/(2*sigma2)) * cstll.pX[j];
        }
    }
    *llr_out = log(tmp0/tmp1);
}

#ifdef INTELMKL
void simulate_awgn(cstll_t cstll, uint64_t *x, double *y, double sigma2, uint64_t nc, VSLStreamStatePtr* rng_stream) {
    double* n = malloc(sizeof(double) * nc);

    vdRngGaussian(VSL_RNG_METHOD_GAUSSIAN_BOXMULLER, rng_stream, nc, n, 0.0, sqrt(sigma2));
    for(size_t i = 0; i < nc; i++) {        
        y[i] = cstll.X[x[i]] + n[i];
    }

    free(n);
}
#endif

double simulate_awgn_legacy(cstll_t cstll, uint64_t *x, double *y, double sigma2, uint64_t nc) {
    double n = 0;
    double Pn = 0;
    double Px = 0;

    for(size_t i = 0; i < nc; i++) {
        n = randn() * sqrt(sigma2);
        Pn += n * n;
        Px += cstll.X[x[i]] * cstll.X[x[i]];
        y[i] = cstll.X[x[i]] + n;
    }

    return Px/Pn;

}

void simulate_bec(uint64_t *x, double* y, double epsilon, uint64_t nc) {

    for(size_t i = 0; i < nc; i++) {
        double val = rand() / (RAND_MAX + 1.0);
        if(val <= epsilon) {
            y[i] = -1.0;
        } else {
            y[i] = (double) x[i];
        }
    }

}


double randn() {
        static double U, V;
        static int phase = 0;
        double Z;

        if(phase == 0) {
                U = (rand() + 1.) / (RAND_MAX + 2.);
                V = rand() / (RAND_MAX + 1.);
                Z = sqrt(-2 * log(U)) * sin(2 * M_PI * V);
        } else
                Z = sqrt(-2 * log(U)) * cos(2 * M_PI * V);

        phase = 1 - phase;

        return Z;
}

int8_t sign(double a) {
    return (a <= 0) ? -1 : 1;
}

double mean(double *x, size_t len) {
    double val = 0;
    for(size_t i = 0; i < len; i++) {
        val += x[i];
    }
    return val/len;
}

uint64_t sum(uint64_t *summands, size_t len) {
    uint64_t val = 0;
    for(size_t i = 0; i < len; i++) {
      val += summands[i];
    }
    return val;
}

double fsum(double *summands, size_t len) {
    double val = 0;
    for(size_t i = 0; i < len; i++) {
      val += summands[i];
    }
    return val;
}

double mean2(double *x, size_t len) {
    double val = 0;
    for(size_t i = 0; i < len; i++) {
        val += x[i] * x[i];
    }
    return val/len;
}

uint8_t check_type(symbols_t* symbols, size_t len, uint64_t* type, size_t type_len) {
    uint64_t *tmp = malloc(sizeof(uint64_t) * type_len);
    for(size_t i = 0; i < type_len; i++) {
        tmp[i] = 0;
    }
    for(size_t i = 0; i < len; i++) {
        tmp[symbols[i]]++;
    }
    for(size_t i = 0; i < type_len; i++) {
        if(type[i] != tmp[i]) {
            free(tmp);
            return 0;
        }
    }
    free(tmp);
    return 1;

}

void reverse_labels(uint16_t* labels, uint16_t* labels_rev, uint16_t M) {
    for(size_t i = 0; i < M; i++) {
        labels_rev[labels[i]] = i;
    }
}

uint64_t ipow(uint64_t base, uint64_t expo) {
    uint64_t res = 1;
    for(size_t i = 0; i < expo; i++) {
        res *= base;
    }
    return res;
}

void shuffle(size_t* idx, size_t len) {
    for (size_t i = 0; i < len - 1; i++) {
        size_t j = i + rand() / (RAND_MAX / (len - i) + 1);
        int t = idx[j];
        idx[j] = idx[i];
        idx[i] = t;
    }
}

void create_nb_interleaver(size_t* pi, size_t permute_len, size_t len) {

    for(size_t i = 0; i < permute_len; i++) {
        pi[i] = i;
    }
    shuffle(pi, permute_len);

    for(size_t i = 0; i < len-permute_len; i++) {
        pi[permute_len+i] = permute_len+i;
        // printf("%d, %d\n", permute_len+i, pi[permute_len+i]);
    }

}

void interleaver(bits_t* input, bits_t* output, size_t* pi, size_t len) {
    for(size_t i = 0; i < len; i++) {
        output[pi[i]] = input[i];
    }
}

void init_double_array(double *arr, size_t len, double val) {
    while(len--) {
        *arr++ = val;
    }
}


/*
 *
 * memory management
 *
 *
 */

void destroy_cstll_t(cstll_t* cstll) {
    free(cstll->X);
    free(cstll->pX);
}

void destroy_dm_t(dm_t* dm) {
    #ifdef SMDM
    free(dm->pA);
    #else
    free(dm->types);
    free(dm->pA);
    #endif
}



/*
 *
 * DEBUGGING functions
 * can be used to inspect data structures
 *
 *
 */

void dec2bin(uint64_t val, uint8_t m) {
    for(size_t i = 0; i < m; i++) {
        printf("%lu", (val>>(m-i-1) & 0x01));
    }
}


void print_cstll_t(cstll_t cstll) {
    double P = 0;
    printf("=========== CSTLL ===========\n");
    printf("M = %d\n", cstll.M);
    for(size_t i = 0; i < cstll.M; i++) {
        printf("%lf: %lf \n", cstll.X[i], cstll.pX[i]);
        P += cstll.X[i] * cstll.X[i] * cstll.pX[i];
    }
    printf("(X.^2.*pX) = %lf\n", P);
    printf("=========== CSTLL: END ===========\n");
}

void print_dm_t(dm_t dm) {
    #ifdef SMDM
    printf("=========== SMDM ===========\n");
    printf("smdm.k: %lu\n", dm.k);
    printf("smdm.n: %lu\n", dm.n);
    printf("smdm.k0: %lu\n", dm.k0);
    printf("smdm.n0: %lu\n", dm.n0);
    printf("smdm.rate: %lf\n", dm.rate);
    printf("=========== END: SMDM ===========\n");
    #else
    printf("=========== CCDM ===========\n");
    printf("ccdm_bits_in: %lu\n", dm.k);
    printf("ccdm_num_types: %lu\n", dm.n);
    printf("ccdm_types: \n");
    int val = 0;
    for(size_t i = 0; i < dm.num_types; i++) {
      printf("\t %lu: %lu\n", i, dm.types[i]);
      val += dm.types[i];
    }
    printf("ccdm_num_symbols_out: %lu\n", dm.n);
    printf("=========== END: CCDM ===========\n");
    #endif
}
