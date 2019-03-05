#include "ldpc_decoder_layered.h"
#include "../function/ldpc_functions.h"


void ldpc_decode_layered(ldpc_code_t* code, double* llr_in, double* l_c2v, double *l_c2v_sum, uint64_t* cn_subset, const uint64_t cn_size)
{
    double* l_v2c;

    size_t* vn;
    size_t* cn;

    size_t vw;
    size_t cw;

    double *f = calloc(code->max_dc, sizeof(double));
    double *b = calloc(code->max_dc, sizeof(double));

    // set to all zero
    l_v2c = calloc(code->nnz, sizeof(double));

    for(size_t i = 0; i < code->nnz; ++i)
        l_c2v[i] = l_c2v_sum[i] - l_c2v[i];

    /* VN node intialization */
    for(size_t i = 0; i < code->nc; i++)
    {
        double tmp = llr_in[i];
        vw = code->vw[i];
        vn = code->vn[i];
        while(vw--)
            tmp += l_c2v[*vn++];

        vn = code->vn[i];
        vw = code->vw[i];
        while(vw--)
        {
            l_v2c[*vn] = tmp - l_c2v[*vn];
            ++vn;
        }
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

    free(l_v2c);

    free(b);
    free(f);
}

uint64_t ldpc_decode_layered_start(ldpc_code_t* code, double *llr_in, double* llr_out, const uint64_t MaxIter)
{
    uint64_t set1[3] = {0,2,4};
    uint64_t set2[3] = {1,3,5};

    size_t* vn;
    size_t vw;

    double lsum[12] = {0.0};
    bits_t bits[6] = {0};

    double *l_c2v_1 = calloc(code->nnz, sizeof (double));
    double *l_c2v_2 = calloc(code->nnz, sizeof (double));


    for (size_t I = 0; I < MaxIter; ++I)
    {

        //parallel
        ldpc_decode_layered(code, llr_in, l_c2v_1, lsum, set1, 3);
        ldpc_decode_layered(code, llr_in, l_c2v_2, lsum, set2, 3);

        //interchange check node messages
        for(int i=0; i<code->nnz; ++i)
        {
            lsum[i] = l_c2v_1[i] + l_c2v_2[i];
            //l_c2v_2[i] = lsum[i] - l_c2v_2[i];
            //l_c2v_1[i] = lsum[i] - l_c2v_1[i];
        }

        // app calculation
        for(size_t i = 0; i < code->nc; i++)
        {
            llr_out[i] = llr_in[i];
            vn = code->vn[i];
            vw = code->vw[i];
            while(vw--)
                llr_out[i] += lsum[*vn++];
        }

        printf("Layerd @Iter: %lu :: \t", I);
        printVectorDouble(llr_out, code->nc);
        for(int i=0; i<code->nc; ++i)
            bits[i] = llr_out[i] <= 0;
        printf("Is codeword: %i \n", is_codeword(*code, bits));


    }


    free(l_c2v_1);free(l_c2v_2);

    return MaxIter;
}
