#include "ldpc_stoppingsets.h"
#include "functions.h"
#include "ldpc_decoder.h"

void lpdc_code_t_stopping_sets(ldpc_code_t *code)
{
    /*
    //for (size_t j = 0; j < code->nc; ++j)
    //for (size_t l = j+1; l < code->nc; ++l)
    size_t j = 0;
    size_t l = 1;
    const size_t I_max = 1; // 10 .. 100
    //const size_t T = 10.0;

    bits_t* bits = calloc(code->nc, sizeof(bits_t));
    double* llr = calloc(code->nc, sizeof(double));
    for (size_t i = 0; i < code->nc; ++i)
        llr[i] = 1.0;

    llr[j] = -10000.0;
    llr[l] = -10000.0;

    for (size_t i = 0; i < I_max; ++i)
    {
        ldpc_decode(*code, llr, &llr, 100, 0, &bits);
    }

    printBits(bits, code->nc);

    free(bits);
    free(llr);
    */
}

void lpdc_code_t_erasure_decoding(ldpc_code_t* code, bits_t* in_bits, bits_t** out_bits)
{
    size_t* epsilon = calloc(1, sizeof(size_t));

    size_t tmp = 0;
    for (size_t i = 0; i < code->nc; ++i)
    {
        if (in_bits[i] > 1) //erasure
        {
            epsilon = realloc(epsilon, (tmp+1) * sizeof(size_t));
            epsilon[tmp] = i;
            ++tmp;
        }
    }

    ldpc_code_t e_code;

    // submatrix of H
    generate_submatrix(code, &e_code, epsilon, tmp);

    print_ldpc_code_t(e_code);

    free(epsilon);
    destroy_ldpc_code_t(&e_code);
}

void generate_submatrix(ldpc_code_t* code, ldpc_code_t* erasure_code, size_t* epsilon, const size_t epsilon_size)
{
    // intialise
    erasure_code->nc = epsilon_size;
    erasure_code->nct = epsilon_size;

    erasure_code->mc = code->mc;
    erasure_code->mct = code->mc;

    // reduce VNs to size of Epsilon
    erasure_code->vw = calloc(epsilon_size, sizeof(size_t));
    for (size_t i = 0; i < epsilon_size; ++i)
        erasure_code->vw[i] = code->vw[epsilon[i]];


    // number of CNs stay the same
    erasure_code->cw = calloc(erasure_code->mc, sizeof(size_t));
    erasure_code->cn = calloc(erasure_code->mc, sizeof(size_t*));
    for (size_t i = 0; i < erasure_code->mc; ++i)
        erasure_code->cn[i] = calloc(0, sizeof(size_t));

    printVector(epsilon, epsilon_size);
    printVector(erasure_code->vw, epsilon_size);

    size_t tmp_rc = 0;
    erasure_code->r = calloc(tmp_rc, sizeof(size_t));
    erasure_code->c = calloc(tmp_rc, sizeof(size_t));

    erasure_code->vn = calloc(epsilon_size, sizeof(size_t*));
    for (size_t i = 0; i < epsilon_size; ++i)
    {
        erasure_code->vn[i] = calloc(erasure_code->vw[i], sizeof(size_t));

        for (size_t j = 0; j < erasure_code->vw[i]; ++j)
        {
            size_t cn_current = code->r[code->vn[epsilon[i]][j]];

            erasure_code->r = realloc(erasure_code->r, (tmp_rc+1) * sizeof(size_t));
            erasure_code->c = realloc(erasure_code->c, (tmp_rc+1) * sizeof(size_t));
            erasure_code->r[tmp_rc] = cn_current;
            erasure_code->c[tmp_rc] = i;

            erasure_code->vn[i][j] = tmp_rc;

            // realloc current CN neighbour & add the current VN
            erasure_code->cn[cn_current] = realloc(erasure_code->cn[cn_current], (erasure_code->cw[cn_current] + 1) * sizeof(size_t));
            erasure_code->cn[cn_current][erasure_code->cw[cn_current]] = tmp_rc;

            // increment weight
            erasure_code->cw[cn_current]++;

            ++tmp_rc;
        }
    }

    erasure_code->nnz = tmp_rc;
}

























