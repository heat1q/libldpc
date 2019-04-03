#include "ldpc_decoder_layered.h"
#include "../function/ldpc_functions.h"


uint64_t ldpc_decode_layered(ldpc_decoder_lyr_t* dec, ldpc_code_t* code, double* llr_in, double* llr_out, const uint64_t MaxIter, const uint8_t early_termination)
{
    size_t* vn;
    size_t* cn;

    size_t vw;
    size_t cw;

    //initialize
    for (size_t i = 0; i < code->nnz; ++i)
    {
        dec->lsum[i] = 0.0;
        for (size_t l = 0; l < code->nl; ++l)
        {
            dec->l_c2v[l][i] = 0.0;
            dec->l_v2c[l][i] = 0.0;
            dec->l_c2v_pre[l][i] = 0.0;
        }
    }

    size_t I = 0;
    while (I < MaxIter)
    {
        for (size_t l = 0; l < code->nl; ++l)
        {
            /* VN node intialization */
            for(size_t i = 0; i < code->nc; i++)
            {
                double tmp = llr_in[i];
                vw = code->vw[i];
                vn = code->vn[i];
                while(vw--)
                    tmp += dec->lsum[*vn++];

                vn = code->vn[i];
                vw = code->vw[i];
                while(vw--)
                {
                    dec->l_v2c[l][*vn] = tmp - dec->l_c2v[l][*vn];
                    ++vn;
                }
            }

            /* CN node processing */
            for(size_t i = 0; i < code->lw[l]; i++)
            {
                cw = code->cw[code->layers[l][i]];
                cn = code->cn[code->layers[l][i]];
                dec->f[l][0] = dec->l_v2c[l][*cn];
                dec->b[l][cw-1] = dec->l_v2c[l][*(cn+cw-1)];
                for(size_t j = 1; j < cw; j++)
                {
                    dec->f[l][j] = jacobian(dec->f[l][j-1], dec->l_v2c[l][*(cn+j)]);
                    dec->b[l][cw-1-j] = jacobian(dec->b[l][cw-j], dec->l_v2c[l][*(cn + cw-j-1)]);
                }

                dec->l_c2v[l][*cn] = dec->b[l][1];
                dec->l_c2v[l][*(cn+cw-1)] = dec->f[l][cw-2];

                for(size_t j = 1; j < cw-1; j++)
                    dec->l_c2v[l][*(cn+j)] = jacobian(dec->f[l][j-1], dec->b[l][j+1]);
            }

            //update the llr sum of layers, by replacing old llr of lyr l with new value
            for (size_t i = 0; i < code->nnz; ++i)
            {
                dec->lsum[i] += dec->l_c2v[l][i] - dec->l_c2v_pre[l][i];
                dec->l_c2v_pre[l][i] = dec->l_c2v[l][i];
            }

            // app calculation
            for(size_t i = 0; i < code->nc; ++i)
            {
                llr_out[i] = llr_in[i];
                vn = code->vn[i];
                vw = code->vw[i];
                while(vw--)
                    llr_out[i] += dec->lsum[*vn++];
                dec->c_out[i] = (llr_out[i] <= 0);
            }

            if (early_termination)
            {
                if (is_codeword(*code, dec->c_out))
                    return I;
            }
        }

        ++I;
    }

    return I;
}

uint8_t layered_dec_setup(ldpc_decoder_lyr_t* dec, ldpc_code_t* code, char* clfile)
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

    //intialize & allocate memory
    dec->lsum = calloc(code->nnz, sizeof(double));
    dec->c_out = calloc(code->nc, sizeof(bits_t));

    dec->l_c2v = calloc(code->nl, sizeof(double*));
    dec->l_v2c = calloc(code->nl, sizeof(double*));
    dec->f = calloc(code->nl, sizeof(double*));
    dec->b = calloc(code->nl, sizeof(double*));
    dec->l_c2v_pre = calloc(code->nl, sizeof(double*));
    for (size_t i = 0; i < code->nl; ++i)
    {
        dec->l_c2v[i] = calloc(code->nnz, sizeof(double));
        dec->l_v2c[i] = calloc(code->nnz, sizeof(double));
        dec->f[i] = calloc(code->max_dc, sizeof(double));
        dec->b[i] = calloc(code->max_dc, sizeof(double));
        dec->l_c2v_pre[i] = calloc(code->nnz, sizeof(double));
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

void destroy_dec(ldpc_decoder_lyr_t *dec, ldpc_code_t* code)
{
    //free memory
    for (size_t i = 0; i < code->nl; ++i)
    {
        free(dec->l_c2v[i]);
        free(dec->l_v2c[i]);
        free(dec->f[i]);
        free(dec->b[i]);
        free(dec->l_c2v_pre[i]);
    }
    free(dec->l_c2v);
    free(dec->l_v2c);
    free(dec->f);
    free(dec->b);
    free(dec->lsum);
    free(dec->l_c2v_pre);
    free(dec->c_out);
}
