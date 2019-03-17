#include "simulation.h"
#include <string.h>
#include <math.h>

using namespace ldpc;

Sim_AWGN_cl::Sim_AWGN_cl(Ldpc_Code_cl *code, const char *simFileName, const char *mapFileName)
{
    setup_sim(code, simFileName, mapFileName);
}

Sim_AWGN_cl::~Sim_AWGN_cl() {}

void Sim_AWGN_cl::setup_sim(Ldpc_Code_cl *code, const char *simFileName, const char *mapFileName)
{
    ldpc_code = code;

    FILE *fp;

    char tmp[1000];
    double tmpd[50];
    double tmpi[1024];
    double m;
    char* tokp;
    size_t i;

    cstll = new cstll_t();

    decoder_terminate_early = 1;

    fp = fopen(simFileName, "r");
    if(!fp)
    {
        printf("can not open sim file\n");
        exit(EXIT_FAILURE);
    }
    fscanf(fp, "name: %256s\n", logfile);
    /*****************************************/
    fscanf(fp, "M: %hu\n", &M);
    cstll->M = M;
    cstll->log2M = static_cast<uint16_t> (log2(cstll->M));
    /*****************************************/
    fscanf(fp, "bits: %hu\n", &bits);
    /*****************************************/
    fscanf(fp, "labels: %[^\n]\n", tmp);
    tokp = strtok(tmp, ", ");
    i = 0;
    while(tokp != nullptr)
    {
        tmpi[i] = atoi(tokp);
        tokp = strtok(nullptr, ", ");
        i++;
    }
    if(i != M)
    {
        printf("error parsing simfile: number of constellation points (%d) does not match label size (%lu)\n", M, i);
        exit(EXIT_FAILURE);
    }

    labels = new uint16_t[M] ();
    labels_rev = new uint16_t[M] ();
    for(size_t j = 0; j < M; j++)
        labels[j] = static_cast<uint16_t> (tmpi[j]);

    //reverse labels
    for(size_t i = 0; i < M; i++)
        labels_rev[labels[i]] = i;

    /*****************************************/
    fscanf(fp, "snrs: %[^\n]\n", tmp);
    tokp = strtok(tmp, ",");
    i = 0;
    while(tokp != nullptr)
    {
        tmpd[i] = static_cast<double> (atof(tokp));
        tokp = strtok(nullptr, ",");
        i++;
    }
    snrs = new double[i];
    for(size_t j = 0; j < i; j++)
        snrs[j] = tmpd[j];
    num_snrs = i;
    /*****************************************/
    fscanf(fp, "max frames: %lu\n", &max_frames);
    /*****************************************/
    fscanf(fp, "min fec: %lu\n", &min_fec);
    /*****************************************/
    fscanf(fp, "bp iter: %lu\n", &bp_iter);

    if(code->nct() % bits != 0)
    {
        printf("Chosen setting m = %hu and channel code with n_c = %lu does not work. Please correct.\n", bits, code->nct());
        exit(EXIT_FAILURE);
    }
    n = code->nct()/bits;

    /* setup constellation */
    cstll->pX = new double[M];
    for(size_t j = 0; j < M; j++) {
        cstll->pX[j] = 1.0/M;
    }
    m = 0;
    cstll->X = new double[M];
    for(size_t j = 0; j < M; j++)
    {
        cstll->X[j] = static_cast<double> (-M+1+2*j);
        m += cstll->X[j] * cstll->X[j] * cstll->pX[j];
    }
    for(size_t j = 0; j < M; j++)
        cstll->X[j] = cstll->X[j]/sqrt(m);

    SE = (static_cast<double>(code->kct())/code->nct()) * bits;

    bit_mapper = new size_t*[bits];
    for(size_t j = 0; j < bits; j++)
        bit_mapper[j] = new size_t[n];

    read_bit_mapping_file(mapFileName);

    bits_pos = new size_t[code->nct()];
    bits_t found_p = 0;
    bits_t found_s = 0;

    size_t idx = 0;
    for(size_t i = 0; i < code->nc(); i++)
    {
        for(size_t j = 0; j < code->num_shorten(); j++)
        {
            if(code->shorten()[j] == i)
            {
                found_s = 1;
                break;
            }
        }
        for(size_t j = 0; j < code->num_puncture(); j++)
        {
            if(code->puncture()[j] == i)
            {
                found_p = 1;
                break;
            }
        }
        if((!found_s) && (!found_p))
        {
            bits_pos[idx] = i;
            idx++;
        }
        else
        {
            found_p = 0;
            found_s = 0;
        }
    }

    fclose(fp);
}

void Sim_AWGN_cl::read_bit_mapping_file(const char* filename)
{
    FILE *fp;

    fp = fopen(filename, "r");

    if(!fp)
    {
        printf("can not open mapping file\n");
        exit(EXIT_FAILURE);
    }

    for(size_t i = 0; i < bits; i++)
    {
        for(size_t j = 0; j < n; j++)
            fscanf(fp, " %lu, ", &(bit_mapper[i][j]));
    }

    fclose(fp);
}

void Sim_AWGN_cl::print_sim()
{
    printf("=========== SIM ===========\n");
    printf("logfile: %s\n", logfile);
    printf("n: %lu\n", n);
    printf("M: %hu\n", M);
    printf("bits: %hu\n", bits);
    printf("SNRs: ");
    for(size_t i = 0; i < num_snrs; i++)
        printf("%lf ", snrs[i]);

    printf("\n");
    printf("labels:\n");
    for(size_t i = 0; i < M; i++)
    {
        printf("\t%lu: ", i);
        dec2bin(labels[i], bits);
        printf("\n");
    }
    printf("max frames: %lu\n", max_frames);
    printf("min fec: %lu\n", min_fec);
    printf("bp iter: %lu\n", bp_iter);
    printf("SE: %.4lf\n", SE);
    // for(size_t i = 0; i < sim.bits; i++) {
    //     printf("Bit Mapping B%lu: ", i);
    //     for(size_t j = 0; j < sim.n; j++) {
    //         printf("%lu, ", sim.bit_mapper[i][j]);
    //     }
    //     printf("\n");
    // }
    printf("=========== SIM: END ===========\n");
}

void Sim_AWGN_cl::destroy_sim()
{
    delete[] snrs;
    delete[] labels;
    delete[] labels_rev;
    for(size_t i = 0; i < bits; i++)
        delete[] bit_mapper[i];

    delete[] bits_pos;
    delete[] bit_mapper;

    //constellation
    delete[] cstll->pX;
    delete[] cstll->X;
    delete cstll;
}

