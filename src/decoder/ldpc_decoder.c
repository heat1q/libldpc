#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../function/scm_functions.h"
#include "ldpc_decoder.h"
#include "../function/ldpc_functions.h"


uint64_t ldpc_decode(ldpc_code_t code, double* llr_in, double** llr_out, uint64_t max_iter, uint8_t early_termination) {
    size_t it;

    double* l_v2c;
    double* l_c2v;

    size_t* vn;
    size_t* cn;

    size_t vw;
    size_t cw;

    bits_t *c_out = calloc(code.nc, sizeof(bits_t));

    // double f[code.max_dc];
    // double b[code.max_dc];
    double *f = malloc(sizeof(double) * code.max_dc);
    double *b = malloc(sizeof(double) * code.max_dc);


    l_v2c = malloc(sizeof(double) * code.nnz);
    l_c2v = malloc(sizeof(double) * code.nnz);

    /* initialize with llrs */
    for(size_t i = 0; i < code.nnz; i++) {
        l_v2c[i] = llr_in[code.c[i]];
    }

    it = 0;
    while(it < max_iter) {
        for(size_t i = 0; i < code.mc; i++) {
            cw = code.cw[i];
            cn = code.cn[i];
            f[0] = l_v2c[*cn];
            b[cw-1] = l_v2c[*(cn+cw-1)];
            for(size_t j = 1; j < cw; j++) {
                f[j] = jacobian(f[j-1], l_v2c[*(cn+j)]);
                b[cw-1-j] = jacobian(b[cw-j], l_v2c[*(cn + cw-j-1)]);
            }

            l_c2v[*cn] = b[1];
            l_c2v[*(cn+cw-1)] = f[cw-2];
            for(size_t j = 1; j < cw-1; j++) {
                l_c2v[*(cn+j)] = jacobian(f[j-1], b[j+1]);
            }
        }

        /* VN node processing */
        for(size_t i = 0; i < code.nc; i++) {
            double tmp = llr_in[i];
            vw = code.vw[i];
            vn = code.vn[i];
            while(vw--) {
                tmp += l_c2v[*vn++];
            }
            vn = code.vn[i];
            vw = code.vw[i];
            while(vw--) {
                l_v2c[*vn] = tmp - l_c2v[*vn];
                vn++;
            }
        }

        // app calculation
        for(size_t i = 0; i < code.nc; i++) {
            (*llr_out)[i] = llr_in[i];
            vn = code.vn[i];
            vw = code.vw[i];
            while(vw--) {
                (*llr_out)[i] += l_c2v[*vn++];
            }
            c_out[i] = ((*llr_out)[i] <= 0);
        }

        it++;

        if(early_termination) {
            if(is_codeword(code, c_out)) {
                break;
            }
        }
    }

    free(l_c2v);
    free(l_v2c);
    free(c_out);

    free(b);
    free(f);

    return it;
}

double jacobian(double L1, double L2) {
#ifdef CN_APPROX_LIN
    return sign(L1) * sign(L2) * fmin(fabs(L1),fabs(L2)) + jacobian_lin_approx(L1+L2) - jacobian_lin_approx(L1-L2);
#elif CN_APPROX_MINSUM
    return sign(L1) * sign(L2) * fmin(fabs(L1), fabs(L2));
#else
    return sign(L1) * sign(L2) * fmin(fabs(L1),fabs(L2)) + log((1+exp(-fabs(L1+L2)))/(1+exp(-fabs(L1-L2))));
#endif
}

double jacobian_lin_approx(double L) {
    double Labs = fabs(L);

    if(Labs < 1.0) {
        return -0.375 * Labs  + 0.6825;
    } else if((Labs >= 1.0) && (Labs < 2.625)) {
        return -0.1875 * Labs + 0.5;
    } else {
        return 0;
    }
}
