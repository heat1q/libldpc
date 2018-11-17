#include "functions.h"

void read_ldpc_file(ldpc_code_t* code, char* filename) {
    FILE *fp;

    fp = fopen(filename, "r");
    if(!fp) {
        printf("can not open codefile %s for reading.\n", filename);
        exit(EXIT_FAILURE);
    }

    fscanf(fp, "nc: %lu\n", &(code->nc));
    fscanf(fp, "mc: %lu\n", &(code->mc));
    fscanf(fp, "nct: %lu\n", &(code->nct));
    fscanf(fp, "mct: %lu\n", &(code->mct));
    fscanf(fp,  "nnz: %lu\n", &(code->nnz));
    code->kc = code->nc-code->mc;
    code->kct = code->nct-code->mct;

    fscanf(fp, "puncture [%lu]: ", &(code->num_puncture));
    code->num_puncture_sys = 0;
    code->num_puncture_par = 0;
    if(code->num_puncture != 0) {
        code->puncture = calloc(code->num_puncture, sizeof(size_t));
        for(size_t i = 0; i < code->num_puncture; i++) {
            fscanf(fp, " %lu ", &(code->puncture[i]));
            if(code->puncture[i] < code->kc) {
                code->num_puncture_sys++;
            } else {
                code->num_puncture_par++;
            }
        }
    } else {
        code->puncture = NULL;
    }
    fscanf(fp, "shorten [%lu]: ", &(code->num_shorten));
    if(code->num_shorten != 0) {
        code->shorten = calloc(code->num_shorten, sizeof(size_t));
        for(size_t i = 0; i < code->num_shorten; i++) {
            fscanf(fp, " %lu ", &(code->shorten[i]));
        }
    } else {
        code->shorten = NULL;
    }


    size_t* cw_tmp;
    size_t* vw_tmp;
    code->cw = calloc(code->mc, sizeof(uint64_t));
    cw_tmp = calloc(code->mc, sizeof(uint64_t));
    code->vw = calloc(code->nc, sizeof(uint64_t));
    vw_tmp = calloc(code->nc, sizeof(uint64_t));
    code->r = calloc(code->nnz, sizeof(size_t));
    code->c = calloc(code->nnz, sizeof(size_t));


    // uint64_t r, c;
    for(size_t i = 0; i < code->nnz; i++) {
        fscanf(fp, "%lu %lu\n", &(code->r[i]), &(code->c[i]));
        // add_ldpc_edge_t(code, i, r, c);
        code->cw[code->r[i]]++;
        code->vw[code->c[i]]++;
    }

    code->cn = calloc(code->mc, sizeof(size_t*));
    for(size_t i = 0; i < code->mc; i++) {
        code->cn[i] = calloc(code->cw[i], sizeof(size_t));
    }
    code->vn = calloc(code->nc, sizeof(size_t*));
    for(size_t i = 0; i < code->nc; i++) {
        code->vn[i] = calloc(code->vw[i], sizeof(size_t));
    }

    for(size_t i = 0; i < code->nnz; i++) {
        code->cn[code->r[i]][cw_tmp[code->r[i]]++] = i;
        code->vn[code->c[i]][vw_tmp[code->c[i]]++] = i;
    }

    free(vw_tmp);
    free(cw_tmp);

    code->max_dc = 0;
    for(size_t i = 0; i < code->mc; i++) {
        if(code->cw[i] > code->max_dc) {
            code->max_dc = code->cw[i];
        }
    }

    code->st_max_size = 0;

#ifdef ENCODE
    code->genmat = calloc(code->kc, sizeof(uint8_t*));
    for(size_t i = 0; i < code->kc; i++) {
        code->genmat[i] = calloc(code->mc, sizeof(uint8_t));
    }
    size_t num;
    size_t tmp;
    for(size_t i = 0; i < code->mc; i++) {
        fscanf(fp, "%lu", &num);
        for(size_t j = 0; j < num; j++) {
            fscanf(fp, "%lu", &tmp);
            code->genmat[tmp][i] = 1;
        }
    }
#endif


    fclose(fp);
}

void print_ldpc_code_t(ldpc_code_t code) {
    printf("=========== LDPC ===========\n");
    printf("nc : %lu\n", code.nc);
    printf("mc : %lu\n", code.mc);
    printf("kc : %lu\n", code.kc);
    printf("nnz : %lu\n", code.nnz);
    printf("nct : %lu\n", code.nct);
    printf("mct : %lu\n", code.mct);
    printf("kct : %lu\n", code.kct);
    printf("num puncture: %lu\n", code.num_puncture);
    printf("num puncture sys: %lu\n", code.num_puncture_sys);
    printf("num puncture par: %lu\n", code.num_puncture_par);
    printf("num shorten: %lu\n", code.num_shorten);
    printf("\n");
    printf("=========== LDPC: END ===========\n");
}


void destroy_ldpc_code_t(ldpc_code_t* code) {
    for(size_t i = 0; i < code->nc; i++) {
        free(code->vn[i]);
    }
    for(size_t i = 0; i < code->mc; i++) {
        free(code->cn[i]);
    }
    free(code->vn);
    free(code->cn);
    free(code->vw);
    free(code->cw);
    free(code->r);
    free(code->c);


    if (code->puncture != NULL)
        free(code->puncture);
    if (code->shorten != NULL)
        free(code->shorten);

    if (code->st_max_size)
    {
        free(code->stw);
        for(size_t i = 0; i < code->st_max_size + 1; i++)
            free(code->st[i]);
        free(code->st);

        free(code->dw);
        for(size_t i = 0; i < code->st_max_size + 1; i++)
            free(code->ds[i]);
        free(code->ds);
    }


#ifdef ENCODE
    for(size_t i = 0; i < code->kc; i++) {
        free(code->genmat[i]);
    }
    free(code->genmat);
#endif
}

void calc_syndrome_c(ldpc_code_t code, bits_t* c, bits_t* s) {
    bits_t b = 0;
    for(size_t i = 0; i < code.mc; i++) {
        b = 0;
        for(size_t j = 0; j < code.cw[i]; j++) {
            b ^= c[code.c[code.cn[i][j]]];
        }
        s[i] = b;
    }
}

bits_t is_codeword(ldpc_code_t code, bits_t* c) {
    bits_t* synd = malloc(sizeof(bits_t) * code.mc);
    bits_t is_codeword = 1;

    calc_syndrome_c(code, c, synd);

    for(size_t i = 0; i < code.mc; i++) {
        if(synd[i] == 1) {
            is_codeword = 0;
            break;
        }
    }

    free(synd);

    return is_codeword;
}

void printVector(void* input, const size_t k)
{
    uint64_t *x = input;
    printf("[");
    for (size_t i = 0; i < k; ++i)
    {
        printf("%lu", x[i]);
        if (i < k - 1)
            printf(" ");
    }
    printf("]\n");
}

void printVectorDouble(double* x, const size_t k)
{
    printf("[");
    for (size_t i = 0; i < k; ++i)
    {
        printf("%.2f", x[i]);
        if (i < k - 1)
            printf(" ");
    }
    printf("]\n");
}

void printBits(bits_t* x, const size_t k)
{
    printf("[");
    for (size_t i = 0; i < k; ++i)
    {
        printf("%i", x[i]);
        if (i < k - 1)
            printf(" ");
    }
    printf("]\n");
}

// just for testing small codes
void generic_codeword_search(ldpc_code_t *code, bits_t** bits, const size_t length, size_t currentbitindex)
{
    for (bits_t b = 0; b < 2; ++b)
    {
        (*bits)[currentbitindex] = b;
        if ((currentbitindex+1) < length)
            generic_codeword_search(code, bits, length, (currentbitindex+1));
        else
        {
            if (is_codeword(*code, *bits))
            {
                printf("====================\n");
                printf("Found valid codeword. \n");
                printBits(*bits, length);
            }
        }
    }
}

size_t NchooseK(const size_t n, const size_t k)
{
    if (k == 0)
        return 1;
    return (n * NchooseK(n - 1, k - 1)) / k;
}
