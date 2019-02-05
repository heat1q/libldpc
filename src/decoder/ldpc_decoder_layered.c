#include "ldpc_decoder_layered.h"
#include "../function/ldpc_functions.h"


uint64_t ldpc_decode_layered(ldpc_code_t* code, double* llr_in, double** llr_out, uint64_t max_iter, uint64_t *vn_subset, const uint64_t vn_size, uint64_t *cn_subset, const uint64_t cn_size)
{
    size_t it;

    double* l_v2c;
    double* l_c2v;

    size_t* vn;
    size_t* cn;

    size_t vw;
    size_t cw;


    // double f[code.max_dc];
    // double b[code.max_dc];
    double *f = malloc(sizeof(double) * code->max_dc);
    double *b = malloc(sizeof(double) * code->max_dc);


    l_v2c = malloc(sizeof(double) * code->nnz);
    l_c2v = malloc(sizeof(double) * code->nnz);

    /* initialize with llrs */
    for(size_t i = 0; i < code->nnz; i++)
        l_v2c[i] = llr_in[code->c[i]];

    it = 0;
    while(it < max_iter)
    {
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

        /* VN node processing */
        for(size_t i = 0; i < vn_size; i++)
        {
            double tmp = llr_in[vn_subset[i]];
            vw = code->vw[vn_subset[i]];
            vn = code->vn[vn_subset[i]];

            while(vw--)
                tmp += l_c2v[*vn++];

            vn = code->vn[vn_subset[i]];
            vw = code->vw[vn_subset[i]];

            while(vw--)
            {
                l_v2c[*vn] = tmp - l_c2v[*vn];
                vn++;
            }
        }

        // app calculation
        for(size_t i = 0; i < vn_size; i++)
        {
            (*llr_out)[vn_subset[i]] = llr_in[vn_subset[i]];
            vn = code->vn[vn_subset[i]];
            vw = code->vw[vn_subset[i]];
            while(vw--)
                (*llr_out)[vn_subset[i]] += l_c2v[*vn++];

        }

        it++;
    }

    free(l_c2v);
    free(l_v2c);

    free(b);
    free(f);

    return it;
}
