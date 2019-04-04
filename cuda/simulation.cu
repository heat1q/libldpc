#include "simulation.cuh"
#include <string.h>
#include <math.h>
#include <exception>
#include <curand.h>
#include <curand_kernel.h>

using namespace std;
using namespace ldpc;

Sim_AWGN_cl::Sim_AWGN_cl(Ldpc_Code_cl* code, const char* simFileName, const char* mapFileName)
: ldpc_code(code)
{
    //init
    cstll = nullptr;
    snrs = nullptr;
    labels = nullptr;
    labels_rev = nullptr;
    bit_mapper = nullptr;
    bits_pos = nullptr;
    ldpc_dec = nullptr;
    
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
        
        decoder_terminate_early = true;
        
        fp = fopen(simFileName, "r");
        if(!fp)
        {
            throw runtime_error("can not open sim file");
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
            labels[j] = tmpi[j];
        
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
        
        //set up decoder
        ldpc_dec = new Ldpc_Decoder_cl(code, bp_iter, decoder_terminate_early, true);
        
        fclose(fp);
    }
    catch(exception &e)
    {
        cout << "Error: " << e.what() << endl;
        destroy();
        
        exit(EXIT_FAILURE);
    }
}

Sim_AWGN_cl::~Sim_AWGN_cl() { destroy(); }


void Sim_AWGN_cl::destroy()
{
    if (snrs != nullptr) { delete[] snrs; }
    if (labels != nullptr) { delete[] labels; }
    if (labels_rev != nullptr) { delete[] labels_rev; }
    if (bits_pos != nullptr) { delete[] bits_pos; }
    if (ldpc_dec != nullptr) { delete ldpc_dec; }
    
    if (bit_mapper != nullptr)
    {
        for(size_t i = 0; i < bits; i++) { delete[] bit_mapper[i]; }
        delete[] bit_mapper;
    }
    
    //constellation
    if (cstll != nullptr)
    {
        if (cstll->pX != nullptr) { delete[] cstll->pX; }
        if (cstll->X != nullptr) { delete[] cstll->X; }
        delete cstll;
    }
}


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

void Sim_AWGN_cl::print()
{
    ldpc_code->print_ldpc_code();
    
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
        if(std::isinf(val) == +1) {
            llrs_out[i] = MAX_LLR;
        } else if(std::isinf(val) == -1) {
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

__host__ __device__ double Sim_AWGN_cl::randn()
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

void Sim_AWGN_cl::encode_all0(uint64_t* x, bits_t* c)
{
    for(size_t i = 0; i < ldpc_code->nct(); i++)
        c[bits_pos[i]] = rand() & 1;
    
    for(size_t i = 0; i < ldpc_code->num_puncture(); i++)
        c[ldpc_code->puncture()[i]] = rand() & 1;
    
    for(size_t i = 0; i < ldpc_code->num_shorten(); i++)
        c[ldpc_code->shorten()[i]] = 0;
    
    
    map_c_to_x(c, x);
}

void Sim_AWGN_cl::map_c_to_x(bits_t* c, size_t* x)
{
    size_t tmp;
    
    for(size_t i = 0; i < n; i++)
    {
        tmp = 0;
        for(size_t j = 0; j < bits; j++)
            tmp += c[bit_mapper[j][i]] << (bits-1-j);
        
        x[i] = labels_rev[tmp];
    }
}


void Sim_AWGN_cl::start()
{
    uint64_t* x;
    double* y;
    bits_t* c;
    double* l_tmp;
    
    double sigma2;
    
    uint64_t frames;
    uint64_t bec = 0;
    uint64_t fec = 0;
    uint64_t iters;
    
    size_t bec_tmp;
    
    FILE* fp = fopen(logfile, "w");
    if(!fp)
    {
        cout << "can not open logfile " << logfile << " for writing" << endl;
        exit(EXIT_FAILURE);
    }
    
    fprintf(fp, "snr fer ber frames avg_iter\n");
    
    
    /*
     * START: SIMULATION PART
     */
    
    for (size_t i = 0; i < num_snrs; ++i)
    {
        bec = 0;
        fec = 0;
        frames = 0;
        iters = 0;
        sigma2 = pow(10, -snrs[i]/10);
        
        /* parallel here */
        {
            x = new uint64_t[n];
            c = new bits_t[ldpc_code->nc()];
            y = new double[n];
            
            do
            {
                encode_all0(x, c);
                
                simulate_awgn(x, y, sigma2);
                
                //puncturing & shortening
                if(ldpc_code->num_puncture() != 0)
                {
                    for(size_t j = 0; j < ldpc_code->num_puncture(); j++) {
                        ldpc_dec->llr_in[ldpc_code->puncture()[j]] = 0;
                    }
                }
                if(ldpc_code->num_shorten() != 0)
                {
                    for(size_t j = 0; j < ldpc_code->num_shorten(); j++) {
                        ldpc_dec->llr_in[ldpc_code->shorten()[j]] = 99999.9;
                    }
                }
                
                double l_tmp[bits];
                for(size_t j = 0; j < n; j++)
                {
                    calc_llrs(y[j], sigma2, l_tmp);
                    for(size_t k = 0; k < bits; k++) {
                        ldpc_dec->llr_in[bit_mapper[k][j]] = l_tmp[k];
                    }
                }
                
                for(size_t j = 0; j < ldpc_code->nc(); j++) {
                    ldpc_dec->llr_in[j] *= (1-2*c[j]);
                }
                
                //decode
                //iters += ldpc_dec->decode_legacy();
                //iters += ldpc_dec->decode_layered_legacy();
                iters += ldpc_dec->decode_layered();
                
                frames++;
                
                bec_tmp = 0;
                for(size_t j = 0; j < ldpc_code->nc(); j++)
                {
                    bec_tmp += (ldpc_dec->llr_out[j] <= 0);
                }
                
                if (bec_tmp > 0)
                {
                    bec += bec_tmp;
                    fec++;
                    
                    printf("FRAME ERROR (%lu/%lu) in frame %lu @SNR = %.3f: BER=%.2e, FER=%.2e\n", fec, min_fec, frames, 10*log10(1/sigma2), (double) bec/(frames*ldpc_code->nc()), (double) fec/frames);
                    
                    //log_error(sim, code, cstll, c, l_out, frames, 10*log10(1/sigma2));
                }
                
            } while (fec < min_fec && frames < max_frames); //end while
            
            delete[] x;
            delete[] c;
            delete[] y;
        }//end parallel
        
        fprintf(fp, "%lf %.3e %.3e %lu %.3e\n", snrs[i], static_cast<double>(fec/frames), static_cast<double>(bec/(frames*ldpc_code->nc())), frames, static_cast<double>(iters/frames));
        fflush(fp);
    }//end for
    
    //legacy decoder
    //delete dec;
    
    fclose(fp);
}

//tmpl fcts need definition in each file?
template<typename T> void ldpc::printVector(T *x, const size_t &l)
{
    cout << "[";
    for (size_t i = 0; i < l-1; ++i)
        cout << x[i] << " ";
    cout << x[l-1] << "]";
}


/*
 *	Cudakernels
 */
__global__ void cudakernel::sim::sim_test(Ldpc_Decoder_cl* dec_mgd)
{
    curandState_t state;
    curand_init(clock64(), 1, 0, &state);
    for (size_t i=0; i<dec_mgd->ldpc_code->nc(); ++i)
    {
        dec_mgd->llr_in[i] = curand_normal(&state);
        dec_mgd->llr_out[i] = 0.0;
    }
}
