#include "ldpc_stoppingsets.h"
#include "functions.h"
#include "ldpc_decoder.h"

#define I_MAX 10 // 10 .. 100
#define THRESH 10
#define REL -10.0

void lpdc_code_t_stopping_sets(ldpc_code_t* code, const size_t MAX_SIZE)
{
    printf("=========== LDPC Stopping Sets ===========\n");
    printf("Maximum Size: %lu \n", MAX_SIZE);

    ldpc_code_t_st_setup(code, MAX_SIZE);

    const size_t t_ges = NchooseK(code->nc, 2);
    size_t time = 0;
    size_t st_found = 0;

    size_t tr_set_size;
    uint8_t tr_stored;
    size_t tmp;

    double* llr;
    bits_t* bits;
    size_t* tr_set;

    char f_w[100]; char f_st[100];
    sprintf(f_w, "st_sets_count.txt");
    sprintf(f_st, "st_sets.txt");
    FILE* file_w = fopen(f_w, "w");
    FILE* file_st = fopen(f_st, "w");

    printf("Result Files: %s %s\n", f_st, f_w);

    struct timespec tstart={0,0}, tend={0,0};
    clock_gettime(CLOCK_MONOTONIC, &tstart);

    #pragma omp parallel for default(none) private(llr, bits, tr_set, tr_set_size, tr_stored, tmp) shared(code, time, file_w, f_w, st_found)
    for (size_t j = 0; j < code->nc; ++j)
    {
        llr = calloc(code->nc, sizeof(double));
        bits = calloc(code->nc, sizeof(bits_t));

        for (size_t l = j+1; l < code->nc; ++l)
        {
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

                    //erasure decoding
                    tr_set_size = lpdc_code_t_erasure_decoding(code, &bits, &tr_set);
                    // when stopping set is found
                    if (tr_set_size && tr_set_size <= MAX_SIZE)
                    {
                        #pragma omp critical
                        {
                            tr_stored = 0;
                            size_t cur_stw = code->stw[tr_set_size];

                            for (size_t a = 0; a < cur_stw; ++a) // for all stored st sets
                            {
                                tmp = 0;
                                for (size_t b = a*tr_set_size; b < (a+1)*tr_set_size; ++b)
                                    tmp += (code->st[tr_set_size][b] == tr_set[b - a*tr_set_size]);
                                if (tmp == tr_set_size)
                                {
                                    tr_stored = 1;
                                    break;
                                }
                            }

                            if (!tr_stored)
                            {
                                //store
                                code->st[tr_set_size] = realloc(code->st[tr_set_size], (cur_stw + 1) * tr_set_size * sizeof(size_t));
                                for (size_t a = 0; a < tr_set_size; ++a)
                                    code->st[tr_set_size][cur_stw * tr_set_size + a] = tr_set[a];

                                ++code->stw[tr_set_size];

                                ++st_found;
                                // distance set
                                /*
                                for (size_t e = 0; e < tr_set_size; ++e)
                                    bits[e] = 0;

                                if (is_codeword(*code, bits))
                                {
                                    code->ds[tr_set_size] = realloc(code->ds[tr_set_size], (code->dw[tr_set_size] + 1) * tr_set_size * sizeof(size_t));
                                    for (size_t a = 0; a < tr_set_size; ++a)
                                        code->ds[tr_set_size][code->dw[tr_set_size] * tr_set_size + a] = tr_set[a];

                                    ++code->dw[tr_set_size];
                                }
                                else
                                {
                                    printf("Not codeword\n");
                                }
                                */
                            }
                            free(tr_set);
                        }
                    }
                    ++t;
                }
            }
            printf("\rProgress: %.2f%% Found: %lu", (double)++time/t_ges *100, st_found);
        }
        free(llr);
        free(bits);
    }
    printf("\n");

    clock_gettime(CLOCK_MONOTONIC, &tend);
    printf("Total Time taken:  %.5f seconds\n", ((double)tend.tv_sec + 1.0e-9*tend.tv_nsec) - ((double)tstart.tv_sec + 1.0e-9*tstart.tv_nsec));

    printVectorToFile(code->stw, MAX_SIZE+1, file_w, 0);

    for (size_t i = 1; i <= MAX_SIZE; ++i)
    {
        fprintf(file_st, "Size=%lu\t# of Stopping Sets=%lu\n", i, code->stw[i]);
        for (size_t j = 0; j < code->stw[i]; ++j)
        {
            printVectorToFile(code->st[i], i, file_st, j*i);
            fprintf(file_st, "\n");
        }
    }

    printf("==========================================\n");
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

    code->stw = calloc(ST_MAX_SIZE + 1, sizeof(size_t));
    code->st = calloc(ST_MAX_SIZE + 1, sizeof(size_t*));
    for (size_t i = 0; i < ST_MAX_SIZE + 1; ++i)
        code->st[i] = calloc(0, sizeof(size_t));

    code->dw = calloc(ST_MAX_SIZE + 1, sizeof(size_t));
    code->ds = calloc(ST_MAX_SIZE + 1, sizeof(size_t*));
    for (size_t i = 0; i < ST_MAX_SIZE + 1; ++i)
        code->ds[i] = calloc(0, sizeof(size_t));
}
