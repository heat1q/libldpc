#define LAYERED_DEC

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <omp.h>

#ifdef INTELMKL
#include <mkl.h>
#endif
#include "function/scm_functions.h"
#include "function/ldpc_functions.h"
#ifdef QUANT
#include "ldpc_decoder_quant.h"
#elif defined LAYERED_DEC
#include "decoder/ldpc_decoder_layered.h"
#else
#include "decoder/ldpc_decoder.h"
#endif



int main(int argc, char* argv[])
{
    /* =========================================================================== */
    uint8_t abort = 0;
    char codeFileName[128];
    char mapFileName[128];
    char simFileName[128];
    char codelayerFileName[128];
    const int numArgs = (argc-1)/2;

    #ifdef LAYERED_DEC
    if (argc == 9)
    {
        for (int i = 0; i < numArgs; ++i)
        {
            if (strcmp(argv[2*i+1], "+code") == 0)
                strcpy(codeFileName, argv[2*i+2]);
            else if (strcmp(argv[2*i+1], "+map") == 0)
                strcpy(mapFileName, argv[2*i+2]);
            else if (strcmp(argv[2*i+1], "+sim") == 0)
                strcpy(simFileName, argv[2*i+2]);
            else if (strcmp(argv[2*i+1], "+cl") == 0)
                strcpy(codelayerFileName, argv[2*i+2]);
            else
                abort = 1;
        }
    }
    else
        abort = 1;

    if (abort)
    {
        printf("======================== LDPC Simulation ========================\n");
        printf("                        (Layered Decoding)                       \n");
        printf("                         Usage Reminder:                         \n");
        printf("         Main +code CodeFile +map CodeMapFile +sim SimFile       \n");
        printf("              +cl CLFile                                         \n");
        printf("                                                                 \n");
        printf("                CodeFile: Name of the code file                  \n");
        printf("               CodeMapFile: Name of mapping file                 \n");
        printf("               SimFile: Name of simulation file                  \n");
        printf("               CLFile: Name of code layer file                   \n");
        printf("=================================================================\n");
        exit(EXIT_FAILURE);
    }
    #else
    if (argc == 7)
    {
        for (int i = 0; i < numArgs; ++i)
        {
            if (strcmp(argv[2*i+1], "+code") == 0)
                strcpy(codeFileName, argv[2*i+2]);
            else if (strcmp(argv[2*i+1], "+map") == 0)
                strcpy(mapFileName, argv[2*i+2]);
            else if (strcmp(argv[2*i+1], "+sim") == 0)
                strcpy(simFileName, argv[2*i+2]);
            else
                abort = 1;
        }
    }
    else
        abort = 1;

    if (abort)
    {
        printf("======================== LDPC Simulation ========================\n");
        printf("                         Usage Reminder:                         \n");
        printf("         Main +code CodeFile +map CodeMapFile +sim SimFile       \n");
        printf("                                                                 \n");
        printf("                CodeFile: Name of the code file                  \n");
        printf("               CodeMapFile: Name of mapping file                 \n");
        printf("               SimFile: Name of simulation file                  \n");
        printf("=================================================================\n");
        exit(EXIT_FAILURE);
    }
    #endif
    /* =========================================================================== */


    ldpc_sim_t sim;
    ldpc_code_t code;
    cstll_t cstll;

    uint64_t* x;
    double* y;
    double* l_in;
    double* l_out;
    bits_t* c;

    double sigma2;

    uint64_t frames;
    uint64_t bec = 0;
    uint64_t fec = 0;
    uint64_t iters;

    FILE *fp = NULL;

    const char* short_opt = "f:i:F:s:t";
    struct option long_opt[] = {
        {"resultfile", required_argument, NULL, 'f'},
        {"iter", required_argument, NULL, 'i'},
        {"numframes", required_argument, NULL, 'F'},
        {"snr", required_argument, NULL, 's'},
        {"noterminate", no_argument, NULL, 't'},
        {NULL, 0, NULL, 0}
    };

    /* read config */

    setup_ldpc_sim(&sim, &code, &cstll, simFileName, codeFileName, mapFileName);

    #ifdef LAYERED_DEC
    ldpc_decoder_lyr_t decoder;
    if(!layered_dec_setup(&decoder, &code, codelayerFileName))
    {
        printf("Can not setup decoder!\n");
        exit(EXIT_FAILURE);
    }
    #endif

    int opt;
    while((opt = getopt_long(argc, argv, short_opt, long_opt, NULL)) != -1) {
        switch(opt) {
            case 'f':
                /* override logfile */
                strncpy(sim.logfile, optarg, MAX_FILENAME_LEN);
                break;
            case 'i':
                /* override number of iterations */
                sim.bp_iter = atoi(optarg);
                break;
            case 'F':
                sim.min_fec = atoi(optarg);
                break;
            case 's':
                /* override snrs, just run this single snr */
                sim.num_snrs = 1;
                sim.snrs[0] = (double) atof(optarg);
                break;
            case 't':
                sim.decoder_terminate_early = 0;
                break;
            default:
                printf("unrecognized parameter\n");
                printf("HELP:\n");
                printf("\t--resultfile (-f): filename of the resultfile\n");
                printf("\t--iter (-i): number of iterations\n");
                printf("\t--snr (-s): only simulate for the given SNR\n");
                printf("\t--numframes (-F): collect at least this number of errors\n");
                printf("\t--noterminate (-t): don't use the syndrome check for early termination\n");
                exit(EXIT_FAILURE);
        }
    }

    print_ldpc_sim_t(sim);
    print_ldpc_code_t(code);
    print_cstll_t(cstll);

    #ifdef INTELMKL
    int num_threads = omp_get_max_threads();
    VSLStreamStatePtr rng_stream[num_threads];
    for(size_t i = 0; i < num_threads; i++) {
        vslNewStream(&rng_stream[i], VSL_BRNG_MT2203+i, 777777);
    }
    #else
    void *rng_stream = NULL;
    #endif

    /*
     * START: SIMULATION PART
     *
     */

    fp = fopen(sim.logfile, "w");
    if(!fp) {
        printf("can not open logfile %s for writing\n", sim.logfile);
        exit(EXIT_FAILURE);
    }
    fprintf(fp, "snr fer ber frames avg_iter\n");

    for(size_t i = 0; i < sim.num_snrs; i++) {
        bec = 0;
        fec = 0;
        frames = 0;
        iters = 0;
        sigma2 = pow(10, -sim.snrs[i]/10);
        #pragma omp parallel default(none) private(c, x, y, l_in, l_out) firstprivate(sim, code, cstll, sigma2) shared(rng_stream, fec, bec, frames) reduction(+:iters)
        {
            #ifdef INTELMKL
            int tid = omp_get_thread_num();
            #endif
           
            x = malloc(sizeof(uint64_t) * sim.n);
            c = malloc(sizeof(bits_t) * code.nc);
            y = malloc(sizeof(double) * sim.n);
            l_in = malloc(sizeof(double) * code.nc);
            l_out = malloc(sizeof(double) * code.nc);
            
            
            while(1) {

                #ifdef ENCODE
                #ifdef INTELMKL
                ldpc_encode(sim, code, x, c, rng_stream[tid]);
                #else
                ldpc_encode_legacy(sim, code, x, c);
                #endif
                #else
                #ifdef INTELMKL
                ldpc_encode_all0(sim, code, x, c, rng_stream[tid]);
                #else
                ldpc_encode_all0_legacy(sim, code, x, c);
                #endif
                #endif

                #ifdef INTELMKL
                simulate_awgn(cstll, x, y, sigma2, sim.n, rng_stream[tid]);
                #else
                simulate_awgn_legacy(cstll, x, y, sigma2, sim.n);
                #endif

                if(code.num_puncture != 0) {
                    for(size_t j = 0; j < code.num_puncture; j++) {
                        l_in[code.puncture[j]] = 0;
                    }
                }
                if(code.num_shorten != 0) {
                    for(size_t j = 0; j < code.num_shorten; j++) {
                        l_in[code.shorten[j]] = 99999.9;
                    }
                }

                double l_tmp[sim.bits];
                for(size_t j = 0; j < sim.n; j++) {
                    calc_llrs(y[j], cstll, sigma2, sim.labels, l_tmp);
                    for(size_t k = 0; k < sim.bits; k++) {
                        l_in[sim.bit_mapper[k][j]] = l_tmp[k];
                    }
                }


                #ifndef ENCODE
                for(size_t j = 0; j < code.nc; j++) {
                    l_in[j] *= (1-2*c[j]);
                }
                #endif

                /* ldpc decoder */
                #if defined QUANT
                iters += ldpc_decode_quant(code, l_in, l_out, sim.bp_iter);
                #elif defined LAYERED_DEC
                iters += ldpc_decode_layered(&decoder, &code, l_in, l_out, sim.bp_iter, sim.decoder_terminate_early);
                #else
                iters += ldpc_decode(code, l_in, l_out, sim.bp_iter, sim.decoder_terminate_early);
                #endif

                // we just processed one frame in this thread, so we increase the
                // global frame count
                #pragma omp atomic update
                frames += 1;

                /* BER and FER tracking */
                size_t bec_tmp = 0;
                #ifdef ENCODE
                for(size_t j = 0; j < code.nc; j++) {
                    bec_tmp += ((l_out[j] <= 0) != c[j]);
                }
                #else
                for(size_t j = 0; j < code.nc; j++)
                {
                    bec_tmp += (l_out[j] <= 0);
                }
                #endif

                if(bec_tmp > 0) {
                    #pragma omp atomic update
                    bec += bec_tmp;
                    #pragma omp atomic update
                    fec += 1;
                    #pragma omp critical 
                    {
                        printf("FRAME ERROR (%lu/%lu) in frame %lu @SNR = %.3f: BER=%.2e, FER=%.2e\n", fec, sim.min_fec, frames, 10*log10(1/sigma2), (double) bec/(frames*code.nc), (double) fec/frames);

                        #ifndef NO_ERR_ANALYSIS
                        log_error(sim, code, cstll, c, l_out, frames, 10*log10(1/sigma2));
                        #endif
                    }
                }

                if(fec >= sim.min_fec || frames >= sim.max_frames) {
                    break;
                }

            } // end: while(1)
            free(c);
            free(x);
            free(y);
            free(l_in);
            free(l_out);
        } // end: pragma omp parallel
        fprintf(fp, "%lf %.3e %.3e %lu %.3e\n", sim.snrs[i], (double) fec/frames, (double) bec/(frames*code.nc), frames, (double) iters / frames);
        fflush(fp);
    }

    /*
     * END: SIMULATION PART
     *
     */

    fclose(fp);

    #ifdef INTELMKL
    for(size_t i = 0; i < num_threads; i++) {
        vslDeleteStream(&rng_stream[i]);
    }
    #endif

    #ifdef LAYERED_DEC
    destroy_dec(&decoder, &code);
    #endif

    destroy_ldpc_code_t(&code);
    destroy_cstll_t(&cstll);
    destroy_ldpc_sim_t(&sim);

    return 0;
}
