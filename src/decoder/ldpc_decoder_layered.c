#include "ldpc_decoder_layered.h"
#include "../function/ldpc_functions.h"


void layered_dec(ldpc_code_t* code, double* llr_in, double* l_c2v, double* l_c2v_sum, double* l_v2c, uint64_t* cn_subset, const uint64_t cn_size, double* f, double* b)
{
    size_t* vn;
    size_t* cn;

    size_t vw;
    size_t cw;

    //for(size_t i = 0; i < code->nnz; ++i)
    //    l_c2v[i] = l_c2v_sum[i] - l_c2v[i];

    /* VN node intialization */
    for(size_t i = 0; i < code->nc; i++)
    {
        double tmp = llr_in[i];
        vw = code->vw[i];
        vn = code->vn[i];
        while(vw--)
        {
            tmp += (l_c2v_sum[*vn] - l_c2v[*vn]);
            vn++;
        }

        vn = code->vn[i];
        vw = code->vw[i];
        while(vw--)
        {
            l_v2c[*vn] = tmp - l_c2v_sum[*vn] + l_c2v[*vn];
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
}

uint64_t ldpc_decode_layered(ldpc_code_t* code, double* llr_in, double* llr_out, const uint64_t MaxIter, const uint8_t early_termination)
{
    size_t* vn;
    size_t vw;

    //intialize & allocate memory
    double* lsum = calloc(code->nnz, sizeof(double));
    bits_t* c_out = calloc(code->nc, sizeof(bits_t));

    double** l_c2v = calloc(code->nl, sizeof(double*));
    double** l_v2c = calloc(code->nl, sizeof(double*));
    double** f = calloc(code->nl, sizeof(double*));
    double** b = calloc(code->nl, sizeof(double*));
    for (size_t i = 0; i < code->nl; ++i)
    {
        l_c2v[i] = calloc(code->nnz, sizeof(double));
        l_v2c[i] = calloc(code->nnz, sizeof(double));
        f[i] = calloc(code->max_dc, sizeof(double));
        b[i] = calloc(code->max_dc, sizeof(double));
    }

    size_t I = 0;
    while (I < MaxIter)
    {
        //parallel
        for (size_t i = 0; i < code->nl; ++i)
            layered_dec(code, llr_in, l_c2v[i], lsum, l_v2c[i], code->layers[i], code->lw[i], f[i], b[i]);

        //interchange check node messages
        for(size_t i = 0; i < code->nnz; ++i)
        {
            lsum[i] = 0.0;
            for (size_t j = 0; j < code->nl; ++j)
                lsum[i] += l_c2v[j][i];
        }

        // app calculation
        for(size_t i = 0; i < code->nc; ++i)
        {
            llr_out[i] = llr_in[i];
            vn = code->vn[i];
            vw = code->vw[i];
            while(vw--)
                llr_out[i] += lsum[*vn++];
            c_out[i] = (llr_out[i] <= 0);
        }

        ++I;

        if (early_termination)
        {
            if (is_codeword(*code, c_out))
                break;
        }
    }

    //free memory
    for (size_t i = 0; i < code->nl; ++i)
    {
        free(l_c2v[i]);
        free(l_v2c[i]);
        free(f[i]);
        free(b[i]);
    }
    free(l_c2v);
    free(l_v2c);
    free(f);
    free(b);
    free(lsum);
    free(c_out);

    return I;
}

uint8_t layered_dec_setup(ldpc_code_t* code, char* clfile)
{
    FILE *fp = fopen(clfile, "r");
    if(!fp)
    {
        printf("Can not open file\n");
        return 0;
    }

    fscanf(fp, "nl: %lu\n", &(code->nl));

    code->lw = calloc(code->nl, sizeof(uint64_t));
    code->layers = calloc(code->nl, sizeof(uint64_t*));

    for (size_t i = 0; i < code->nl; ++i)
    {
        fscanf(fp, "cn[i]: %lu\n", &(code->lw[i]));
        code->layers[i] = calloc(code->lw[i], sizeof(uint64_t));
        for (size_t j = 0; j < code->lw[i]; ++j)
            fscanf(fp, "%lu\n", &(code->layers[i][j]));
    }

    printf("=========== LDPC LAYERS ===========\n");
    printf("nl: %lu\n", code->nl);
    for (size_t i = 0; i < code->nl; ++i)
    {
        printf("cn[%lu]: %lu\n", i, code->lw[i]);
        printVector(code->layers[i], code->lw[i]);
    }
    printf("========= LDPC LAYERS: END ========\n");

    fclose(fp);

    return 1;
}
