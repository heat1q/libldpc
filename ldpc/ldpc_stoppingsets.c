#include "ldpc_stoppingsets.h"
#include "functions.h"
#include "ldpc_decoder.h"


#define I_MAX 10 // 10 .. 100
#define THRESH 10
#define REL -10.0

void lpdc_code_t_stopping_sets(ldpc_code_t *code)
{
    for (size_t j = 0; j < code->nc; ++j)
    {
        for (size_t l = j+1; l < code->nc; ++l)
        {
            double* llr = calloc(code->nc, sizeof(double));
            bits_t* bits = calloc(code->nc, sizeof(bits_t));
            for (size_t i = 0; i < code->nc; ++i)
                llr[i] = 1.0;

            llr[j] = REL; llr[l] = REL;

            for (size_t i = 0; i < I_MAX; ++i)
            {
                // run the belief propagation decoder for one interation
                ldpc_decode(*code, llr, &llr, 1, 0);

                // set the threshold
                int64_t t = -1;

                while (t <= THRESH)
                {
                    // declare all bits for which R_h <= t as erasure & all others as 0
                    for (size_t h = 0; h < code->nc; ++h)
                        bits[h] = (uint8_t)((llr[h] <= t)*2);

                    size_t *tr_set;
                    //erasure decoding
                    const size_t tr_set_size = lpdc_code_t_erasure_decoding(code, &bits, &tr_set);
                    // when stopping set is found

                    if (tr_set_size && tr_set_size <= code->st_max_size)
                    {
                        uint8_t tr_stored = 0;
                        for (size_t a = 0; a < code->stw[tr_set_size-2]; ++a) // for all stored st sets
                        {
                            size_t tmp = 0;
                            for (size_t b = a*tr_set_size; b < (a+1)*tr_set_size; ++b)
                                tmp += (code->st[tr_set_size-2][b] == tr_set[b - a*tr_set_size]);
                            if (tmp == tr_set_size)
                            {
                                tr_stored = 1;
                                break;
                            }
                        }
                        if (!tr_stored)
                        {
                            //store
                            code->st[tr_set_size-2] = realloc(code->st[tr_set_size-2], (code->stw[tr_set_size-2] + 1) * tr_set_size * sizeof(size_t));
                            for (size_t a = 0; a < tr_set_size; ++a)
                                code->st[tr_set_size-2][code->stw[tr_set_size-2] * tr_set_size + a] = tr_set[a];

                            ++code->stw[tr_set_size-2];

                            // stopping set found
                            printf("Stopping set found at VN nodes: ");
                            printVector(tr_set, tr_set_size);
                            //printBits(bits, code->nc);

                            //printf("IsCodeword: %u \n", is_codeword(*code, bits));
                        }

                        free(tr_set);
                    }
                    ++t;
                }
            }

            free(llr);
            free(bits);
        }
    }
}

size_t lpdc_code_t_erasure_decoding(ldpc_code_t* code, bits_t** in_bits, size_t** set)
{
    while (1)
    {
        size_t tmp_e = 0;
        size_t* epsilon = calloc(tmp_e, sizeof(size_t));

        for (size_t i = 0; i < code->nc; ++i)
        {
            if ((*in_bits)[i] > 1) //erasure
            {
                epsilon = realloc(epsilon, (tmp_e+1) * sizeof(size_t));
                epsilon[tmp_e] = i;
                ++tmp_e;
            }
        }

        if (!tmp_e) // all erasures removed
        {
            free(epsilon);
            *set = NULL;
            return 0;
        }

        ldpc_code_t e_code;

        // submatrix of H
        generate_submatrix(code, &e_code, epsilon, tmp_e);

        size_t tmp_r = 0;
        size_t* rows = calloc(tmp_r, sizeof(size_t)); // set of rows with weight 1
        // find the set of rows in H_e with single 1-component
        for (size_t j = 0; j < e_code.mc; ++j)
        {
            if (e_code.cw[j] == 1)
            {
                rows = realloc(rows, (tmp_r+2) * sizeof(size_t));
                rows[tmp_r] = j; // denotes the row of H
                rows[tmp_r+1] = epsilon[e_code.c[e_code.cn[j][0]]]; // gives the position of erased bit

                tmp_r += 2;
            }
        }

        if (!tmp_r) // stop; erasure cannot be recoverd
        {
            *set = epsilon;
            //free(epsilon);
            free(rows);
            destroy_ldpc_code_t(&e_code);
            return tmp_e;
        }

        for (size_t t = 0; t < tmp_r; t+=2)
        {
            // parity-checks
            bits_t bit = 0;
            for (size_t n = 0; n < code->cw[rows[t]]; ++n)
            {
                size_t vn_index = code->c[code->cn[rows[t]][n]];

                bit ^= (*in_bits)[vn_index];
            }

            //recover bit
            (*in_bits)[rows[t+1]] = bit ^ (*in_bits)[rows[t+1]];
        }

        free(epsilon);
        free(rows);
        destroy_ldpc_code_t(&e_code);
    }
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

    erasure_code->puncture = NULL;
    erasure_code->shorten = NULL;
    erasure_code->st_max_size = 0;
}

void ldpc_code_t_st_setup(ldpc_code_t* code, const size_t ST_MAX_SIZE)
{
    code->st_max_size = ST_MAX_SIZE;

    code->stw = calloc(ST_MAX_SIZE - 1, sizeof(size_t));
    code->st = calloc(ST_MAX_SIZE - 1, sizeof(size_t*));
    for (size_t i = 0; i < ST_MAX_SIZE - 1; ++i)
        code->st[i] = calloc(0, sizeof(size_t));
}
