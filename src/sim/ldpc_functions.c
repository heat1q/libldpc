#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "scm_functions.h"
#include "ldpc_functions.h"
#ifdef PAS
#ifdef SMDM
#include "smdmlib_MB.h"
#else
#include "ccdmlib.h"
#endif
#endif
#ifdef INTELMKL
#include <mkl.h>
#endif

void encode(ldpc_code_t code, bits_t* c) {

    bits_t bit;
    for(size_t i = 0; i < code.mc; i++) {
        bit = 0;
        for(size_t j = 0; j < code.kc; j++) {
            bit = bit ^ c[j] * code.genmat[j][i];
        }
        c[code.kc+i] = bit;
    }
}

void calc_syndrome_llr(ldpc_code_t code, double* llr_out, bits_t* s) {
    bits_t b = 0;
    for(size_t i = 0; i < code.mc; i++) {
        b = 0;
        for(size_t j = 0; j < code.cw[i]; j++) {
            b ^= (llr_out[code.c[code.cn[i][j]]] <= 0);
        }
        s[i] = b;
    }
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

size_t analyze_codeword(ldpc_code_t code, bits_t* c, size_t* one_pos) {    
    size_t weight = 0;
    size_t k = 0;

    for(size_t i = 0; i < code.nc; i++) {
        weight += c[i];
        if(c[i] == 1) {
            one_pos[k] = i;
            k++;
        }
    }

    return weight;
}

void map_c_to_x(ldpc_sim_t sim, bits_t* c, size_t* x) {
    size_t tmp;

    for(size_t i = 0; i < sim.n; i++) {
        tmp = 0;
        for(size_t j = 0; j < sim.bits; j++) {
            tmp += c[sim.bit_mapper[j][i]] << (sim.bits-1-j);
        }
        x[i] = sim.labels_rev[tmp];
    }
}

#ifdef INTELMKL
void ldpc_encode_all0(ldpc_sim_t sim, ldpc_code_t code, uint64_t* x, bits_t* c, VSLStreamStatePtr* rng_stream) {
    uint64_t tmp;

    int* bits = malloc(sizeof(int) * (code.nct + code.num_puncture));
    viRngBernoulli(VSL_RNG_METHOD_BERNOULLI_ICDF, rng_stream, (code.nct + code.num_puncture), bits, 0.5);

    for(size_t i = 0; i < code.nct; i++) {
        c[sim.bits_pos[i]] = (bits_t) bits[i];
    }
    for(size_t i = 0; i < code.num_puncture; i++) {
        c[code.puncture[i]] =  (bits_t) bits[code.nct + i];
    }
    for(size_t i = 0; i < code.num_shorten; i++) {
        c[code.shorten[i]] = 0;
    }

    map_c_to_x(sim, c, x);

    free(bits);
}
#endif

void ldpc_encode_all0_legacy(ldpc_sim_t sim, ldpc_code_t code, uint64_t* x, bits_t* c) {

    for(size_t i = 0; i < code.nct; i++) {
        c[sim.bits_pos[i]] = rand() & 1;
    }
    for(size_t i = 0; i < code.num_puncture; i++) {
        c[code.puncture[i]] = rand() & 1;
    }
    for(size_t i = 0; i < code.num_shorten; i++) {
        c[code.shorten[i]] = 0;
    }

    map_c_to_x(sim, c, x);

}

#ifdef INTELMKL
void ldpc_encode(ldpc_sim_t sim, ldpc_code_t code, uint64_t* x, bits_t* c, VSLStreamStatePtr* rng_stream) {
    uint64_t tmp;

    int* bits = malloc(sizeof(int) * code.kc);
    viRngBernoulli(VSL_RNG_METHOD_BERNOULLI_ICDF, rng_stream, code.kc, bits, 0.5);

    for(size_t i = 0; i < code.kc; i++) {
        c[i] = (bits_t) bits[i];
    }
    for(size_t i = 0; i < code.num_shorten; i++) {
        c[code.shorten[i]] = 0;
    }

    encode(code, c);

    map_c_to_x(sim, c, x);

    free(bits);
}
#endif

void ldpc_encode_legacy(ldpc_sim_t sim, ldpc_code_t code, uint64_t* x, bits_t* c) {
    uint64_t tmp;


    for(size_t i = 0; i < code.kc; i++) {
        c[i] = rand() & 1;
    }
    for(size_t i = 0; i < code.num_shorten; i++) {
        c[code.shorten[i]] = 0;
    }

    encode(code, c);

    map_c_to_x(sim, c, x);
}

#ifdef PAS
void ldpc_pas_encode_legacy(ldpc_sim_t sim, ldpc_code_t code, dm_t dm, uint64_t *x, bits_t* c, bits_t* data, uint8_t* symbols) {
    uint64_t tmp;

    for(size_t i = 0; i < dm.k; i++) {
        data[i] = rand() & 1;
    }

#ifdef SMDM
    vecEncode(data, dm.k, symbols, dm.n, dm.k0, dm.n0, 2);
#else
    encodeCCDM(data, dm.k, symbols, dm.types, dm.num_types, 31);
#endif

    for(size_t i = 0; i < sim.n; i++) {
        for(size_t j = 0; j < sim.bits-1; j++) {
            c[sim.bit_mapper[j+1][i]] = (sim.labels[sim.M/2 + symbols[i]] >> (sim.bits - 2 - j)) & 1;
        }
    }

    for(size_t i = 0; i < code.num_shorten; i++) {
        c[code.shorten[i]] = 0;
    }

    for(size_t i = 0; i < code.num_puncture && i < code.kc; i++) {
        c[code.puncture[i]] = rand() & 1;
    }

    // add the gamma bits
    uint64_t num_gamma_bits = code.kct - sim.n * (sim.bits-1);
    for(size_t i = 0; i < num_gamma_bits; i++) {
        c[sim.bit_mapper[0][i]] = rand() & 1;
    }

    encode(code, c);

    map_c_to_x(sim, c, x);

}

#ifdef INTELMKL
void ldpc_pas_encode(ldpc_sim_t sim, ldpc_code_t code, dm_t dm, uint64_t *x, bits_t* c, bits_t* data, uint8_t* symbols, VSLStreamStatePtr* rng_stream) {
    uint64_t tmp;
    size_t num_gamma_bits = code.kct - sim.n * (sim.bits-1);
    size_t offset = 0;

    int* data1 = malloc(sizeof(int) * (dm.k + code.num_puncture_sys + num_gamma_bits));
    viRngBernoulli(VSL_RNG_METHOD_BERNOULLI_ICDF, rng_stream, dm.k + code.num_puncture_sys + num_gamma_bits, data1, 0.5);

    for(size_t i = 0; i < dm.k; i++) {
        data[i] = (bits_t) data1[i];
    }
    offset = dm.k;

#ifdef SMDM
    vecEncode(data, dm.k, symbols, dm.n, dm.k0, dm.n0, 2);
#else
    encodeCCDM(data, dm.k, symbols, dm.types, dm.num_types, 31);
#endif

    for(size_t i = 0; i < sim.n; i++) {
        for(size_t j = 0; j < sim.bits-1; j++) {
            c[sim.bit_mapper[j+1][i]] = (sim.labels[sim.M/2 + symbols[i]] >> (sim.bits - 2 - j)) & 1;
        }
    }

    for(size_t i = 0; i < code.num_shorten; i++) {
        c[code.shorten[i]] = 0;
    }

    // this deals with punctured bits in the systematic part
    for(size_t i = 0; i < code.num_puncture_sys; i++) {
        c[code.puncture[i]] = data1[offset+i];
    }
    offset += code.num_puncture_sys;

    // add the gamma bits
    viRngBernoulli(VSL_RNG_METHOD_BERNOULLI_ICDF, rng_stream, num_gamma_bits, data1, 0.5);
    for(size_t i = 0; i < num_gamma_bits; i++) {
        c[sim.bit_mapper[0][i]] = data1[offset + i];
    }

    encode(code, c);

    map_c_to_x(sim, c, x);

    free(data1);
}
#endif

void ldpc_pas_encode_all0_legacy(ldpc_sim_t sim, ldpc_code_t code, dm_t dm, uint64_t *x, bits_t* c, bits_t* data, uint8_t* symbols) {
    uint64_t tmp;

    for(size_t i = 0; i < dm.k; i++) {
        data[i] = rand() & 1;
    }

#ifdef SMDM
    vecEncode(data, dm.k, symbols, dm.n, dm.k0, dm.n0, 2);
#else
    encodeCCDM(data, dm.k, symbols, dm.types, dm.num_types, 31);
#endif

    for(size_t i = 0; i < sim.n; i++) {
        for(size_t j = 0; j < sim.bits-1; j++) {
            c[sim.bit_mapper[j+1][i]] = (sim.labels[sim.M/2 + symbols[i]] >> (sim.bits - 2 - j)) & 1;
        }
    }

    for(size_t i = 0; i < code.num_shorten; i++) {
        c[code.shorten[i]] = 0;
    }

    for(size_t i = 0; i < code.num_puncture && i < code.kc; i++) {
        c[code.puncture[i]] = rand() & 1;
    }

    // fake the sign bits
    for(size_t i = 0; i < sim.n; i++) {
        c[sim.bit_mapper[0][i]] = rand() & 1;
    }

    map_c_to_x(sim, c, x);
}

#ifdef INTELMKL
void ldpc_pas_encode_all0(ldpc_sim_t sim, ldpc_code_t code, dm_t dm, uint64_t *x, bits_t* c, bits_t* data, uint8_t* symbols, VSLStreamStatePtr* rng_stream) {

    size_t offset = 0;

    int *data1 = malloc(sizeof(int) * (dm.k + code.num_puncture + sim.n));
    viRngBernoulli(VSL_RNG_METHOD_BERNOULLI_ICDF, rng_stream, dm.k + code.num_puncture + sim.n, data1, 0.5);

    for(size_t i = 0; i < dm.k; i++) {
        data[i] = (bits_t) data1[i];
    }
    offset = dm.k;

#ifdef SMDM
    vecEncode(data, dm.k, symbols, dm.n, dm.k0, dm.n0, 2);
#else
    encodeCCDM(data, dm.k, symbols, dm.types, dm.num_types, 31);
#endif

    for(size_t i = 0; i < sim.n; i++) {
        for(size_t j = 0; j < sim.bits-1; j++) {
            c[sim.bit_mapper[j+1][i]] = (sim.labels[sim.M/2 + symbols[i]] >> (sim.bits - 2 - j)) & 1;
        }
    }

    for(size_t i = 0; i < code.num_shorten; i++) {
        c[code.shorten[i]] = 0;
    }

    for(size_t i = 0; i < code.num_puncture; i++) {
        c[code.puncture[i]] = data1[offset + i];
    }
    offset += code.num_puncture;

    // fake the sign bits
    for(size_t i = 0; i < sim.n; i++) {
        c[sim.bit_mapper[0][i]] = data1[offset + i];
    }

    map_c_to_x(sim, c, x);

    free(data1);
}
#endif
#endif

void setup_ldpc_sim_pas(ldpc_sim_t *sim, ldpc_code_t* code, cstll_t* cstll, dm_t *dm, char *simfile, char *ldpcfile, char* mappingfile) {
    FILE *fp;

    char tmp[1000];
    double tmpd[50];
    double tmpi[1024];
    double m;
    char* tokp;
    size_t i;

    read_ldpc_file(code, ldpcfile);

    sim->decoder_terminate_early = 1;

    fp = fopen(simfile, "r");
    if(!fp) {
        printf("can not open sim file\n");
        exit(EXIT_FAILURE);
    }
    fscanf(fp, "name: %256s\n", sim->logfile);
    /*****************************************/
    fscanf(fp, "M: %hu\n", &(sim->M));
    cstll->M = sim->M;
    cstll->log2M = log2(cstll->M);
    /*****************************************/
    fscanf(fp, "bits: %hu\n", &(sim->bits));
    if((code->nct % sim->bits) != 0) {
        printf("Chosen setting Bits = %hu and channel code with n_c = %lu does not work. Please correct.\n", sim->bits, code->nct);
        exit(EXIT_FAILURE);
    }
    sim->n = code->nct/sim->bits;
    /*****************************************/
    fscanf(fp, "labels: %[^\n]\n", tmp);
    tokp = strtok(tmp, ", ");
    i = 0;
    while(tokp != NULL) {
        tmpi[i] = atoi(tokp);
        tokp = strtok(NULL, ", ");
        i++;
    }
    if(i != sim->M) {
        printf("error parsing simfile: number of constellation points (%d) does not match label size (%lu)\n", sim->M, i);
        exit(EXIT_FAILURE);
    }
    sim->labels = malloc(sizeof(uint16_t) * sim->M);
    sim->labels_rev = malloc(sizeof(uint16_t) * sim->M);
    for(size_t j = 0; j < sim->M; j++) {
        sim->labels[j] = tmpi[j];
    }
    reverse_labels(&(sim->labels[0]), &(sim->labels_rev[0]), sim->M);

    /*****************************************/
    
    /* setup DM info */
    /*****************************************/
    dm->n = sim->n;
    
    #ifdef SMDM
    fscanf(fp, "smdm_k: %lu\n", &(dm->k0));
    fscanf(fp, "smdm_n: %lu\n", &(dm->n0));
    
    if((sim->n % dm->n0) != 0) {
        printf("Number of channel uses %lu must be a multiple of the matcher output length %lu\n.", sim->n, dm->n0);
        exit(EXIT_FAILURE);
    }

    dm->rate= ((double) dm->k0/dm->n0);
    dm->k =  (size_t) (dm->rate * dm->n);
    
    fscanf(fp, "smdm_pA: %[^\n]\n", tmp);
    tokp = strtok(tmp, ", ");
    i = 0;
    while(tokp != NULL) {
        tmpd[i] =  atof(tokp);
        tokp = strtok(NULL, ", ");
        i++;
    }
    dm->pA = malloc(sizeof(double) * i);
    for(size_t j = 0; j < sim->M/2; j++) {
        dm->pA[j] = tmpd[j];
    }
    #else
    fscanf(fp, "ccdm_bits_in: %lu\n", &(dm->k));
    fscanf(fp, "ccdm_types: %[^\n]\n", tmp);
    tokp = strtok(tmp, ", ");
    i = 0;
    while(tokp != NULL) {
        tmpi[i] =  atoi(tokp);
        tokp = strtok(NULL, ", ");
        i++;
    }
    if(i != sim->M/2) {
        printf("error parsing simfile: number of constellation points (%hu) does not match half the number of types (%lu)\n", sim->M, i);
        exit(EXIT_FAILURE);
    }
    dm->num_types = i;
    dm->types = malloc(sizeof(uint64_t) * dm->num_types);
    dm->pA = malloc(sizeof(double) * dm->num_types);
    for(size_t j = 0; j < sim->M/2; j++) {
        dm->types[j] = tmpi[j];
        dm->pA[j] = (double) dm->types[j]/dm->n;
    }
    dm->n = sum(dm->types, dm->num_types);
    #endif
    /*****************************************/
    fscanf(fp, "snrs: %[^\n]\n", tmp);
    tokp = strtok(tmp, ",");
    i = 0;
    while(tokp != NULL) {
        tmpd[i] = (double) atof(tokp);
        tokp = strtok(NULL, ",");
        i++;
    }
    sim->snrs = malloc(sizeof(double) * i);
    for(size_t j = 0; j < i; j++) {
        sim->snrs[j] = tmpd[j];
    }
    sim->num_snrs = i;
    /*****************************************/
    fscanf(fp, "max frames: %lu\n", &(sim->max_frames));
    /*****************************************/
    fscanf(fp, "min fec: %lu\n", &(sim->min_fec));
    /*****************************************/
    fscanf(fp, "bp iter: %lu\n", &(sim->bp_iter));

    /* setup constellation */
    cstll->pX = malloc(sizeof(double) * sim->M);
    for(size_t j = 0; j < sim->M; j++) {
        cstll->pX[j] = dm->pA[(j < sim->M/2) ? (sim->M/2-1-(j%(sim->M/2))) : (j % (sim->M/2))];
        cstll->pX[j] /= 2;
    }
    m = 0;
    cstll->X = malloc(sizeof(double)*sim->M);
    for(size_t j = 0; j < sim->M; j++) {
        cstll->X[j] = (double) -sim->M+1+2*j;
        m += cstll->X[j] * cstll->X[j] * cstll->pX[j];
    }
    for(size_t j = 0; j < sim->M; j++) {
        cstll->X[j] = cstll->X[j]/sqrt(m);
    }

    sim->SE = (double) dm->k/dm->n + 1 - (1 - (double) code->kct/code->nct) * sim->bits;

    sim->bit_mapper = malloc(sizeof(size_t*) * sim->bits);
    for(size_t j = 0; j < sim->bits; j++) {
        sim->bit_mapper[j] = malloc(sizeof(size_t) * sim->n);
    }
    read_bit_mapping_file(mappingfile, sim->bit_mapper, sim->bits, sim->n);

    sim->bits_pos = malloc(sizeof(size_t) * (code->nct));
    bits_t found_p = 0;
    bits_t found_s = 0;

    size_t idx = 0;
    for(size_t i = 0; i < code->nc; i++) {
      // printf("%lu\n", idx);
        for(size_t j = 0; j < code->num_shorten; j++) {
            if(code->shorten[j] == i) {
                found_s = 1;
                break;
            }
        }
        for(size_t j = 0; j < code->num_puncture; j++) {
            if(code->puncture[j] == i) {
                found_p = 1;
                break;
            }
        }
        if((!found_s) && (!found_p)) {
            sim->bits_pos[idx] = i;
            idx++;
        } else {
            found_p = 0;
            found_s = 0;
        }
    }

    fclose(fp);
}

void setup_ldpc_sim(ldpc_sim_t *sim, ldpc_code_t* code, cstll_t* cstll, char *simfile, char *ldpcfile, char* mappingfile) {
    FILE *fp;

    char tmp[1000];
    double tmpd[50];
    double tmpi[1024];
    double m;
    char* tokp;
    size_t i;

    read_ldpc_file(code, ldpcfile);

    sim->decoder_terminate_early = 1;

    fp = fopen(simfile, "r");
    if(!fp) {
        printf("can not open sim file\n");
        exit(EXIT_FAILURE);
    }
    fscanf(fp, "name: %256s\n", sim->logfile);
    /*****************************************/
    fscanf(fp, "M: %hu\n", &(sim->M));
    cstll->M = sim->M;
    cstll->log2M = log2(cstll->M);
    /*****************************************/
    fscanf(fp, "bits: %hu\n", &(sim->bits));
    /*****************************************/
    fscanf(fp, "labels: %[^\n]\n", tmp);
    tokp = strtok(tmp, ", ");
    i = 0;
    while(tokp != NULL) {
        tmpi[i] = atoi(tokp);
        tokp = strtok(NULL, ", ");
        i++;
    }
    if(i != sim->M) {
        printf("error parsing simfile: number of constellation points (%d) does not match label size (%lu)\n", sim->M, i);
        exit(EXIT_FAILURE);
    }
    sim->labels = malloc(sizeof(uint16_t) * sim->M);
    sim->labels_rev = malloc(sizeof(uint16_t) * sim->M);
    for(size_t j = 0; j < sim->M; j++) {
        sim->labels[j] = tmpi[j];
    }
    reverse_labels(sim->labels, sim->labels_rev, sim->M);
    /*****************************************/
    fscanf(fp, "snrs: %[^\n]\n", tmp);
    tokp = strtok(tmp, ",");
    i = 0;
    while(tokp != NULL) {
        tmpd[i] = (double) atof(tokp);
        tokp = strtok(NULL, ",");
        i++;
    }
    sim->snrs = malloc(sizeof(double) * i);
    for(size_t j = 0; j < i; j++) {
        sim->snrs[j] = tmpd[j];
    }
    sim->num_snrs = i;
    /*****************************************/
    fscanf(fp, "max frames: %lu\n", &(sim->max_frames));
    /*****************************************/
    fscanf(fp, "min fec: %lu\n", &(sim->min_fec));
    /*****************************************/
    fscanf(fp, "bp iter: %lu\n", &(sim->bp_iter));

    if(code->nct % sim->bits != 0) {
        printf("Chosen setting m = %hu and channel code with n_c = %lu does not work. Please correct.\n", sim->bits, code->nct);
        exit(EXIT_FAILURE);
    }
    sim->n = code->nct/sim->bits;

    /* setup constellation */
    cstll->pX = malloc(sizeof(double) * sim->M);
    for(size_t j = 0; j < sim->M; j++) {
        cstll->pX[j] = 1.0/sim->M;
    }
    m = 0;
    cstll->X = malloc(sizeof(double)*sim->M);
    for(size_t j = 0; j < sim->M; j++) {
        cstll->X[j] = (double) -sim->M+1+2*j;
        m += cstll->X[j] * cstll->X[j] * cstll->pX[j];
    }
    for(size_t j = 0; j < sim->M; j++) {
        cstll->X[j] = cstll->X[j]/sqrt(m);
    }

    sim->SE = (((double) code->kct)/code->nct) * sim->bits;

    sim->bit_mapper = malloc(sizeof(size_t*) * sim->bits);
    for(size_t j = 0; j < sim->bits; j++) {
        sim->bit_mapper[j] = malloc(sizeof(size_t) * sim->n);
    }
    
    read_bit_mapping_file(mappingfile, sim->bit_mapper, sim->bits, sim->n);

    sim->bits_pos = malloc(sizeof(size_t) * (code->nct));
    bits_t found_p = 0;
    bits_t found_s = 0;

    size_t idx = 0;
    for(size_t i = 0; i < code->nc; i++) {
      // printf("%lu\n", idx);
        for(size_t j = 0; j < code->num_shorten; j++) {
            if(code->shorten[j] == i) {
                found_s = 1;
                break;
            }
        }
        for(size_t j = 0; j < code->num_puncture; j++) {
            if(code->puncture[j] == i) {
                found_p = 1;
                break;
            }
        }
        if((!found_s) && (!found_p)) {
            sim->bits_pos[idx] = i;
            idx++;
        } else {
            found_p = 0;
            found_s = 0;
        }
    }

    fclose(fp);
}

int read_ldpc_file(ldpc_code_t* code, const char* filename) {
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
        code->puncture = malloc(sizeof(size_t) * code->num_puncture);
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
        code->shorten = malloc(sizeof(size_t) * code->num_shorten);
        for(size_t i = 0; i < code->num_shorten; i++) {
            fscanf(fp, " %lu ", &(code->shorten[i]));
        }
    } else {
        code->shorten = NULL;
    }

    size_t* cw_tmp;
    size_t* vw_tmp;
    code->cw = malloc(sizeof(uint64_t) * code->mc);
    cw_tmp = malloc(sizeof(uint64_t) * code->mc);
    code->vw = malloc(sizeof(uint64_t) * code->nc);
    vw_tmp = malloc(sizeof(uint64_t) * code->nc);
    code->r = malloc(sizeof(size_t) * code->nnz);
    code->c = malloc(sizeof(size_t) * code->nnz);

    for(size_t i = 0; i < code->mc; i++) {
        code->cw[i] = 0;
        cw_tmp[i] = 0;
    }
    for(size_t i = 0; i < code->nc; i++) {
        code->vw[i] = 0;
        vw_tmp[i] = 0;
    }

    // uint64_t r, c;
    for(size_t i = 0; i < code->nnz; i++) {
        fscanf(fp, "%lu %lu\n", &(code->r[i]), &(code->c[i]));
        // add_ldpc_edge_t(code, i, r, c);
        code->cw[code->r[i]]++;
        code->vw[code->c[i]]++;
    }

    code->cn = malloc(sizeof(size_t*) * code->mc);
    for(size_t i = 0; i < code->mc; i++) {
        code->cn[i] = malloc(sizeof(size_t) * code->cw[i]);
    }
    code->vn = malloc(sizeof(size_t*) * code->nc);
    for(size_t i = 0; i < code->nc; i++) {
        code->vn[i] = malloc(sizeof(size_t) * code->vw[i]);
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


#ifdef ENCODE
    code->genmat = malloc(sizeof(uint8_t*) * code->kc);
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
    return 0;
}

void read_bit_mapping_file(char* filename, size_t** bit_mapping, size_t bits, size_t n) {
    FILE *fp;

    fp = fopen(filename, "r");

    if(!fp) {
      printf("can not open mapping file\n");
      exit(EXIT_FAILURE);
    }

    for(size_t i = 0; i < bits; i++) {
        for(size_t j = 0; j < n; j++) {
            fscanf(fp, " %lu, ", &(bit_mapping[i][j]));
        }
    }

    fclose(fp);
}

void log_error(ldpc_sim_t sim, ldpc_code_t code, cstll_t cstll, bits_t *c, double *lout, size_t frame_num, double snr) {
    char errors_file[MAX_FILENAME_LEN];
    snprintf(errors_file, MAX_FILENAME_LEN, "errors_%s", sim.logfile);
    
    FILE *fp = fopen(errors_file, "a+");
    if(!fp) {
        printf("can not open error log file.\n");
        exit(EXIT_FAILURE);
    }
    
    /* calculation of syndrome and failed syndrome checks */
    bits_t *synd = malloc(sizeof(bits_t) * code.mc);
    calc_syndrome_llr(code, lout, synd);
    size_t synd_weight = 0;
    for(size_t i = 0; i < code.mc; i++) {
        synd_weight += (size_t) synd[i];
    }
    size_t* failed_checks_idx = malloc(sizeof(size_t) * synd_weight);
    size_t j = 0;
    for(size_t i = 0; i < code.mc; i++) {
        if(synd[i] == 1) {
            failed_checks_idx[j++] = i;
        }
    }

    /* calculation of failed codeword bits */
    size_t cw_dis = 0;
    for(size_t i = 0; i < code.nc; i++) {
        #ifdef ENCODE
        cw_dis += ((lout[i] <= 0) != c[i]);
        #else
        cw_dis += ((lout[i] <= 0) != 0);
        #endif
    }
    size_t* x = malloc(sizeof(size_t) * sim.n);
    size_t* xhat = malloc(sizeof(size_t) * sim.n);
    bits_t* chat = malloc(sizeof(bits_t) * code.nc);
    for(size_t i = 0; i < code.nc; i++) {
        chat[i] = (bits_t) (lout[i] <= 0);
    }
    map_c_to_x(sim, c, x);
    map_c_to_x(sim, chat, xhat);
    double cw_dis_euc = 0;
    for(size_t i = 0; i < sim.n; i++) {
        #ifdef ENCODE
        cw_dis_euc += (cstll.X[x[i]] - cstll.X[xhat[i]]) * (cstll.X[x[i]] - cstll.X[xhat[i]]);
        #else
        cw_dis_euc += (cstll.X[0] - cstll.X[xhat[i]]) * (cstll.X[0] - cstll.X[xhat[i]]);
        #endif
    }
    size_t* failed_bits_idx = malloc(sizeof(size_t) * cw_dis);
    j = 0;
    for(size_t i = 0; i < code.nc; i++) {
        #ifdef ENCODE
        if(chat[i] != c[i]) {
            failed_bits_idx[j++] = i;
        }
        #else
        if(chat[i] != 0) {
            failed_bits_idx[j++] = i;
        }
        #endif 
    }

    /* print results in file */
    fprintf(fp, "SNR: %.2f -- frame: %lu -- is codeword: %d -- dE(c,chat): %.3f -- dH(c,chat): %lu | ", snr, frame_num, synd_weight == 0, cw_dis_euc, cw_dis);
    for(size_t i = 0; i < cw_dis; i++) {
        fprintf(fp, "%lu ", failed_bits_idx[i]);
    }
    fprintf(fp, " -- ");
    fprintf(fp, "synd weight: %lu | ", synd_weight);
    for(size_t i = 0; i < synd_weight; i++) {
        fprintf(fp, "%lu ", failed_checks_idx[i]);
    }
    fprintf(fp, "\n");
    fclose(fp);

    /* cleanup */
    free(synd);
    free(failed_checks_idx);
    free(failed_bits_idx);
    free(chat);
    free(xhat);
    free(x);

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
    free(code->puncture);
    free(code->shorten);
#ifdef ENCODE
    for(size_t i = 0; i < code->kc; i++) {
      free(code->genmat[i]);
    }
    free(code->genmat);
#endif
}

void destroy_ldpc_sim_t(ldpc_sim_t *sim) {
    free(sim->snrs);
    free(sim->labels);
    free(sim->labels_rev);
    for(size_t i = 0; i < sim->bits; i++) {
        free(sim->bit_mapper[i]);
    }
    free(sim->bits_pos);
    free(sim->bit_mapper);
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
    printf("max dc : %lu\n", code.max_dc);
    printf("num puncture: %lu\n", code.num_puncture);
    printf("num puncture sys: %lu\n", code.num_puncture_sys);
    printf("num puncture par: %lu\n", code.num_puncture_par);
    printf("num shorten: %lu\n", code.num_shorten);
    printf("\n");
    printf("=========== LDPC: END ===========\n");
}


void print_ldpc_sim_t(ldpc_sim_t sim) {
    printf("=========== SIM ===========\n");
    printf("logfile: %s\n", sim.logfile);
    printf("n: %lu\n", sim.n);
    printf("M: %hu\n", sim.M);
    printf("bits: %hu\n", sim.bits);
    printf("SNRs: ");
    for(size_t i = 0; i < sim.num_snrs; i++) {
        printf("%lf ", sim.snrs[i]);
    }
    printf("\n");
    printf("labels:\n");
    for(size_t i = 0; i < sim.M; i++) {
        printf("\t%lu: ", i);
        dec2bin(sim.labels[i], sim.bits);
        printf("\n");
    }
    printf("max frames: %lu\n", sim.max_frames);
    printf("min fec: %lu\n", sim.min_fec);
    printf("bp iter: %lu\n", sim.bp_iter);
    printf("SE: %.4lf\n", sim.SE);
    // for(size_t i = 0; i < sim.bits; i++) {
    //     printf("Bit Mapping B%lu: ", i);
    //     for(size_t j = 0; j < sim.n; j++) {
    //         printf("%lu, ", sim.bit_mapper[i][j]);
    //     }
    //     printf("\n");
    // }
    printf("=========== SIM: END ===========\n");
}
