#include "ldpc_cycles.h"

void cycle_ldpc_code_t(ldpc_code_t* code)
{
    printf("Counting cycles of LDPC code... ");

    const uint64_t edge_count = code->nnz;
    const size_t length = code->nnz * 2*(code->girth);
    size_t ex_msg_cache = 1; // default message size

    uint64_t* counter = calloc(code->girth, sizeof(uint64_t));
    uint64_t** edge_msg = calloc(length, sizeof(uint64_t*)); // allocate the array for all messages

    for (size_t i = 0; i < length; ++i)
        edge_msg[i] = calloc(ex_msg_cache, sizeof(uint64_t)); // allocate all messages

    for (size_t index = 0; index < code->nc; ++index) //initialise at all nodes u_i
    {
        size_t k = 1;

        const size_t ex_msg_size = code->vw[index]; // degree of initialisation node u_i

        if (ex_msg_cache < ex_msg_size) // reallocation
        {
            ex_msg_cache = ex_msg_size; // cache now equals the msg size

            for (size_t i = 0; i < length; ++i)
            {
                edge_msg[i] = realloc(edge_msg[i], ex_msg_cache * sizeof(uint64_t));
                for (size_t j = 0; j < ex_msg_cache; ++j) // set all entries to zero
                    edge_msg[i][j] = 0;
            }
        }
        else
        {
            for (size_t i = 0; i < length; ++i) // set all entries to zero
            {
                for (size_t j = 0; j < ex_msg_size; ++j)
                    edge_msg[i][j] = 0;
            }
        }

        //Initialisation
        for (size_t i = 0; i < ex_msg_size; ++i) // for all edges on node u_index
        {
            const size_t edge = code->vn[index][i];

            //set data array
            edge_msg[edge + (2*k-2)*edge_count][i] = 1;

            //printf("t=%lu (k=1) :: Message at edge e%lu: ", (2*k-2), edge);
            //printVector(edge_msg[edge + (2*k-2)*edge_count], ex_msg_size);
        }

        //printf("\n");

        for (; k < code->girth; ++k)
        {
            //Message Passing from W
            for (size_t j = 0; j < code->mc; ++j) // forall nodes wj
            {
                for (size_t a = 0; a < code->cw[j]; ++a) // forall neighbors of w_j
                {
                    if (code->c[code->cn[j][a]] >= index) // dont send to inactive nodes
                    {
                        for (size_t b = 0; b < code->cw[j]; ++b) // multiply messages
                        {
                            if (b != a)
                                ex_msg_t_add(edge_msg[code->cn[j][a] + (2*k-1)*edge_count], edge_msg[code->cn[j][b] + (2*k-2)*edge_count], ex_msg_size);
                        }
                        //printf("t=%lu (k=%lu) :: Message at edge e%lu: ", 2*k-1, k, code->cn[j][a]);
                        //printVector(edge_msg[code->cn[j][a] + (2*k-1)*edge_count], ex_msg_size);
                    }
                }
            }

            //printf("\n");

            // Counting Cycles
            uint64_t local = 0;
            for (size_t l = 0; l < ex_msg_size; ++l)
                ex_msg_t_sum(&local, edge_msg[code->vn[index][l] + (2*k-1)*edge_count], ex_msg_size, l);

            if (k < code->girth-1)
            {
                // Message Passing from U
                for (size_t i = index+1; i < code->nc; ++i) // disable inactive nodes
                {
                    for (size_t a = 0; a < code->vw[i]; ++a)
                    {
                        for (size_t b = 0; b < code->vw[i]; ++b)
                        {
                            if (b != a)
                                ex_msg_t_add(edge_msg[code->vn[i][a] + 2*k*edge_count], edge_msg[code->vn[i][b] + (2*k-1)*edge_count], ex_msg_size);
                        }
                        //printf("t=%lu (k=%lu) :: Message at edge e%lu: ", 2*k, k+1, code->vn[i][a]);
                        //printVector(edge_msg[code->vn[i][a] + 2*k*edge_count], ex_msg_size);
                    }
                }
                //printf("\n");
            }

            counter[k] += local/2;
        }

        printf("\rCounting cycles of LDPC code...  %.2f%%", (double)index/code->nc *100);
    }

    printf("\rCounting cycles of LDPC code...  100%% Completed.");
    printf("\n");
    printf("=========== LDPC Cycles ===========\n");

    //print to file
    FILE *fp;
    char fname[100];
    sprintf(fname, "cycles_%lux%lu.txt", code->nc, code->mc);
    fp = fopen(fname, "w");
    fprintf(fp, "length count\n");
    printf("length count\n");
    for (size_t i = 2; i < code->girth; ++i)
    {
        fprintf(fp, "%lu %lu\n", 2*i, counter[i]);
        printf("%lu %lu\n", 2*i, counter[i]);
    }

    printf("===================================\n");
    printf("Results saved in %s.\n", fname);

    for (size_t i = 0; i < length; ++i)
        free(edge_msg[i]);

    free(edge_msg);
}

// calculate girth of ldpc code
void girth_ldpc_code_t(ldpc_code_t* code)
{
    printf("Calculating girth of LDPC code... ");

    uint64_t girth = 0;

    const uint64_t edge_count = code->nnz;
    const size_t length = code->nnz * 2;

    size_t ex_var_size = 1;
    size_t ex_msg_cache = 1;

    uint64_t** edge_msg = calloc(length*ex_var_size, sizeof(uint64_t*));

    for (size_t i = 0; i < length*ex_var_size; ++i)
        edge_msg[i] = calloc(ex_msg_cache, sizeof(uint64_t));

    for (size_t index = 0; index < code->nc; ++index)
    {
        size_t k = 1;
        const size_t ex_msg_size = code->vw[index];

        if (ex_msg_cache < ex_msg_size) // reallocation
        {
            ex_msg_cache = ex_msg_size; // cache now equals the msg size

            for (size_t i = 0; i < length*ex_var_size; ++i)
            {
                edge_msg[i] = realloc(edge_msg[i], ex_msg_cache * sizeof(uint64_t));
                for (size_t j = 0; j < ex_msg_cache; ++j) // set all entries to zero
                    edge_msg[i][j] = 0;
            }
        }
        else
        {
            for (size_t i = 0; i < length*ex_var_size; ++i) // set all entries to zero
            {
                for (size_t j = 0; j < ex_msg_size; ++j)
                    edge_msg[i][j] = 0;
            }
        }

        //Initialisation
        for (size_t i = 0; i < ex_msg_size; ++i) // for all edges on node u_index
            edge_msg[code->vn[index][i] + (2*k-2)*edge_count][i] = 1;

        while (1)
        {
            if (2*k > ex_var_size) // reallocation
            {
                ++ex_var_size;
                edge_msg = realloc(edge_msg, length*ex_var_size * sizeof(uint64_t*));

                for (size_t i = length*(ex_var_size-1); i < length*ex_var_size; ++i)
                    edge_msg[i] = calloc(ex_msg_size, sizeof(uint64_t));
            }

            //Message Passing from W
            for (size_t j = 0; j < code->mc; ++j)
            {
                for (size_t a = 0; a < code->cw[j]; ++a)
                {
                    if (code->c[code->cn[j][a]] >= index)
                    {
                        for (size_t b = 0; b < code->cw[j]; ++b)
                        {
                            if (b != a)
                                ex_msg_t_add(edge_msg[code->cn[j][a] + (2*k-1)*edge_count], edge_msg[code->cn[j][b] + (2*k-2)*edge_count], ex_msg_size);
                        }
                    }
                }
            }

            uint64_t local = 0;
            for (size_t l = 0; l < ex_msg_size; ++l)
                ex_msg_t_sum(&local, edge_msg[code->vn[index][l] + (2*k-1)*edge_count], ex_msg_size, l);

            if (local > 0 || k > 10) // dont get caught in infinte loop
            {
                if (girth > 2*k || girth == 0)
                    girth = 2*k;

                break;
            }

            // Message Passing from U
            for (size_t i = index+1; i < code->nc; ++i)
            {
                for (size_t a = 0; a < code->vw[i]; ++a)
                {
                    for (size_t b = 0; b < code->vw[i]; ++b)
                    {
                        if (b != a)
                            ex_msg_t_add(edge_msg[code->vn[i][a] + 2*k*edge_count], edge_msg[code->vn[i][b] + (2*k-1)*edge_count], ex_msg_size);
                    }
                }
            }

            ++k;

            if (2*k >= girth && girth > 0)
                break;
        }

        printf("\rCalculating girth of LDPC code...(local=%lu)  %.2f%%", girth, (double)index/code->nc *100);
    }

    for (size_t i = 0; i < length*ex_var_size; ++i)
        free(edge_msg[i]);

    free(edge_msg);

    printf("\rCalculating girth of LDPC code...(girth=%lu)  100%% Completed.", girth);
    printf("\n");

    code->girth = girth;
}

void printVector(uint64_t* x, const size_t k)
{
    printf("[");
    for (size_t i = 0; i < k; ++i)
    {
        printf("%lu", x[i]);
        if (i < k - 1)
            printf(" ");
    }
    printf("]\n");
}

void ex_msg_t_add(uint64_t* result, uint64_t* x, const size_t length)
{
    for (size_t i = 0; i < length; ++i)
        result[i] += x[i];
}

void ex_msg_t_sum(uint64_t* result, uint64_t* data, const size_t length, const size_t j)
{
    for (size_t i = 0; i < length; ++i)
        *result += data[i];

    *result -= data[j]; // dont count cycles with tails
}
