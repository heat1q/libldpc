#include "simulation.h"
#include <string.h>
#include <math.h>
#include <exception>

using namespace std;
using namespace ldpc;

Sim_AWGN_cl::Sim_AWGN_cl(Ldpc_Code_cl *code, const char *simFileName, const char *mapFileName)
{
    //init
    ldpc_code = code;
    cstll = nullptr;
    snrs = nullptr;
    labels = nullptr;
    labels_rev = nullptr;
    bit_mapper = nullptr;
    bits_pos = nullptr;

    try
    {
        FILE *fp;

        char tmp[1000];
        double tmpd[50];
        double tmpi[1024];
        double m;
        char* tokp;
        size_t i;

        cstll = new cstll_t();
        cstll->pX = nullptr;
        cstll->X = nullptr;

        decoder_terminate_early = 1;

        fp = fopen(simFileName, "r");
        if(!fp)
        {
            throw runtime_error("can not open sim file");
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
            throw runtime_error("error parsing simfile: number of constellation points does not match label size");
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
            throw runtime_error("Chosen setting m with n_c does not work. Please correct.");
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
    catch(exception &e)
    {
        cout << "Error: " << e.what() << endl;

        destroy_sim();
    }
}

Sim_AWGN_cl::~Sim_AWGN_cl() { destroy_sim(); }



void Sim_AWGN_cl::read_bit_mapping_file(const char* filename)
{
    FILE *fp;

    fp = fopen(filename, "r");

    if(!fp)
    {
        throw runtime_error("can not open mapping file");
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
    if (snrs != nullptr)
        delete[] snrs;

    if (labels != nullptr)
        delete[] labels;

    if (labels_rev != nullptr)
        delete[] labels_rev;

    if (bit_mapper != nullptr)
    {
        for(size_t i = 0; i < bits; i++)
            delete[] bit_mapper[i];
        delete[] bit_mapper;
    }

    if (bits_pos != nullptr)
        delete[] bits_pos;

    //constellation
    if (cstll != nullptr)
    {
        if (cstll->pX != nullptr)
            delete[] cstll->pX;
        if (cstll->X != nullptr)
            delete[] cstll->X;
        delete cstll;
    }
}

void Sim_AWGN_cl::calc_llrs(const double &y, const double &sigma2, double *llrs_out)
{
    double tmp0, tmp1;

    for(size_t i = 0; i < cstll->log2M; i++)
    {
        tmp0 = 0.0;
        tmp1 = 0.0;
        for(size_t j = 0; j < cstll->M; j++) {
            if(labels[j] & (1 << (cstll->log2M-1-i))) {
                tmp1 += exp(-(y-cstll->X[j])*(y-cstll->X[j])/(2*sigma2)) * cstll->pX[j];
            } else {
                tmp0 += exp(-(y-cstll->X[j])*(y-cstll->X[j])/(2*sigma2)) * cstll->pX[j];
            }
        }
        double val = log(tmp0/tmp1);
        // check usually required when PAS is used with large constellations
        // and severely shaped distributions
        if(isinf(val) == +1) {
            llrs_out[i] = MAX_LLR;
        } else if(isinf(val) == -1) {
            llrs_out[i] = MIN_LLR;
        } else {
            llrs_out[i] = val;
        }
    }
}

double Sim_AWGN_cl::simulate_awgn(uint64_t *x, double *y, const double &sigma2)
{
    double a = 0;
    double Pn = 0;
    double Px = 0;

    for (size_t i = 0; i < n; i++)
    {
        a = randn() * sqrt(sigma2);
        Pn += a * a;
        Px += cstll->X[x[i]] * cstll->X[x[i]];
        y[i] = cstll->X[x[i]] + a;
    }

    return Px/Pn;
}

double Sim_AWGN_cl::randn()
{
    static double U, V;
    static int phase = 0;
    double Z;

    if (phase == 0)
    {
        U = (rand() + 1.) / (RAND_MAX + 2.);
        V = rand() / (RAND_MAX + 1.);
        Z = sqrt(-2 * log(U)) * sin(2 * M_PI * V);
    }
    else
        Z = sqrt(-2 * log(U)) * cos(2 * M_PI * V);

    phase = 1 - phase;

    return Z;
}
