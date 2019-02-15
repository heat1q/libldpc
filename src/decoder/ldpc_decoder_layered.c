#include "ldpc_decoder_layered.h"
#include "../function/ldpc_functions.h"


void ldpc_decode_layered(ldpc_code_t* code, double** llr, uint64_t *cn_subset, const uint64_t cn_size, const size_t SubIter)
{
    double* l_v2c;
    double* l_c2v;

    size_t* vn;
    size_t* cn;

    size_t vw;
    size_t cw;

    double *f = calloc(code->max_dc, sizeof(double));
    double *b = calloc(code->max_dc, sizeof(double));

    // set to all zero
    l_v2c = calloc(code->nnz, sizeof(double));
    l_c2v = calloc(code->nnz, sizeof(double));

    for (size_t it = 0; it < SubIter; ++it)
    {

        /* VN node intialization */
        for(size_t i = 0; i < code->nnz; i++) //TODO - optimize
        {
            l_v2c[i] = (*llr)[code->c[i]];
            l_c2v[i] = 0.0;
        }

        /* CN node processing */
        for(size_t i = 0; i < cn_size; i++)
        {
            cw = code->cw[cn_subset[i]];
            cn = code->cn[cn_subset[i]];
            f[0] = l_v2c[*cn];
            b[cw-1] = l_v2c[*(cn+cw-1)];
            for(size_t j = 1; j < cw; j++)
            {
                f[j] = jacobian(f[j-1], l_v2c[*(cn+j)]);
                b[cw-1-j] = jacobian(b[cw-j], l_v2c[*(cn + cw-j-1)]);
            }

            l_c2v[*cn] = b[1];
            l_c2v[*(cn+cw-1)] = f[cw-2];

            for(size_t j = 1; j < cw-1; j++)
                l_c2v[*(cn+j)] = jacobian(f[j-1], b[j+1]);
        }

        // app calculation
        for(size_t i = 0; i < code->nc; ++i)
        {
            vn = code->vn[i];
            vw = code->vw[i];
            while(vw--)
                (*llr)[i] += l_c2v[*vn++];
        }
    }

    free(l_c2v);
    free(l_v2c);

    free(b);
    free(f);
}

uint64_t ldpc_decode_layered_init(ldpc_code_t* code, double** llr, const uint64_t MaxIter)
{
    for (uint64_t i = 0; i < MaxIter; ++i)
    {
        // for test code
        ldpc_decode_layered(code, llr, LUT_testcode_CN_subset_1, 400, 1);
        ldpc_decode_layered(code, llr, LUT_testcode_CN_subset_2, 400, 1);
        ldpc_decode_layered(code, llr, LUT_testcode_CN_subset_3, 400, 1);
    }

    return MaxIter;
}
