#include "simulation.h"
#include <string.h>

Sim::Sim()
{

}

Sim::Sim(LDPC *code, const char *simFileName, const char *mapFileName)
{
    setup_sim(code, simFileName, mapFileName);
}

bool Sim::setup_sim(LDPC *code, const char *simFileName, const char *mapFileName)
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
        return false;
    }
    fscanf(fp, "name: %256s\n", logfile);
    /*****************************************/
    fscanf(fp, "M: %hu\n", &M);
    cstll->M = M;
    cstll->log2M = log2(cstll->M);
    /*****************************************/
    fscanf(fp, "bits: %hu\n", &bits);
    /*****************************************/
    fscanf(fp, "labels: %[^\n]\n", tmp);
    tokp = strtok(tmp, ", ");
    i = 0;
    while(tokp != NULL)
    {
        tmpi[i] = atoi(tokp);
        tokp = strtok(NULL, ", ");
        i++;
    }
    if(i != M)
    {
        printf("error parsing simfile: number of constellation points (%d) does not match label size (%lu)\n", M, i);
        return false;
    }

    labels = new uint16_t[M] ();
    labels_rev = new uint16_t[M] ();
    for(size_t j = 0; j < M; j++)
        labels[j] = tmpi[j];

    //reverse labels
    for(size_t i = 0; i < M; i++)
        labels_rev[labels[i]] = i;

    /*****************************************/
    fscanf(fp, "snrs: %[^\n]\n", tmp);
    tokp = strtok(tmp, ",");
    i = 0;
    while(tokp != NULL)
    {
        tmpd[i] = (double) atof(tokp);
        tokp = strtok(NULL, ",");
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
        return false;
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
        cstll->X[j] = (double) -M+1+2*j;
        m += cstll->X[j] * cstll->X[j] * cstll->pX[j];
    }
    for(size_t j = 0; j < M; j++)
        cstll->X[j] = cstll->X[j]/sqrt(m);

    SE = (((double) code->kct())/code->nct()) * bits;

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

    return true;
}

void Sim::read_bit_mapping_file(const char* filename)
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

void Sim::destroy_sim()
{
    delete[] snrs;
    delete[] labels;
    delete[] labels_rev;
    for(size_t i = 0; i < sim->bits; i++)
        delete[] bit_mapper[i];

    delete[] bits_pos;
    delete[] bit_mapper;

    //constellation
    delete[] cstll->pX;
    delete[] cstll->X;
    delete cstll;
}
