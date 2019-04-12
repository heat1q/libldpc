#include "ldpcsim.h"

#include <math.h>
#include <curand.h>
#include <curand_kernel.h>

using namespace ldpc;

ldpc_sim::ldpc_sim(ldpc_code* code, const char* simFileName, const char* mapFileName)
: mLdpcCode(code)
{
    //init
    cstll = nullptr;
    snrs = nullptr;
    labels = nullptr;
    labels_rev = nullptr;
    bit_mapper = nullptr;
    bits_pos = nullptr;
    mLdpcDecoder = nullptr;

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
            throw std::runtime_error("can not open sim file");
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
            throw std::runtime_error("error parsing simfile: number of constellation points does not match label size");
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
            throw std::runtime_error("Chosen setting m with n_c does not work. Please correct.");
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
        mLdpcDecoder = new ldpc_decoder(code, bp_iter, decoder_terminate_early);
        fclose(fp);
    }
    catch(std::exception &e)
    {
        std::cout << "Error: " << e.what() << "\n";
        destroy();

        exit(EXIT_FAILURE);
    }
}

ldpc_sim::~ldpc_sim() { destroy(); }


void ldpc_sim::destroy()
{
    if (snrs != nullptr) { delete[] snrs; }
    if (labels != nullptr) { delete[] labels; }
    if (labels_rev != nullptr) { delete[] labels_rev; }
    if (bits_pos != nullptr) { delete[] bits_pos; }
    if (mLdpcDecoder != nullptr) { delete mLdpcDecoder; }

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


void ldpc_sim::read_bit_mapping_file(const char* filename)
{
    FILE *fp;

    fp = fopen(filename, "r");

    if(!fp)
    {
        throw std::runtime_error("can not open mapping file");
    }

    for(size_t i = 0; i < bits; i++)
    {
        for(size_t j = 0; j < n; j++)
            fscanf(fp, " %lu, ", &(bit_mapper[i][j]));
    }

    fclose(fp);
}

void ldpc_sim::print()
{
    mLdpcCode->print();

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


void ldpc_sim::calc_llrs(const double &y, const double &sigma2, double *llrs_out)
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

double ldpc_sim::simulate_awgn(uint64_t *x, double *y, const double &sigma2)
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

double ldpc_sim::randn()
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

void ldpc_sim::encode_all0(uint64_t* x, bits_t* c)
{
    for(size_t i = 0; i < mLdpcCode->nct(); i++)
        c[bits_pos[i]] = rand() & 1;

    for(size_t i = 0; i < mLdpcCode->num_puncture(); i++)
        c[mLdpcCode->puncture()[i]] = rand() & 1;

    for(size_t i = 0; i < mLdpcCode->num_shorten(); i++)
        c[mLdpcCode->shorten()[i]] = 0;


    map_c_to_x(c, x);
}

void ldpc_sim::map_c_to_x(bits_t* c, size_t* x)
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


void ldpc_sim::start()
{
    uint64_t* x = new uint64_t[n];
    double* y = new double[n];
    bits_t* c = new bits_t[mLdpcCode->nc()];
    double* l_tmp = new double[bits];

    double sigma2;
    uint64_t frames;
    uint64_t bec = 0;
    uint64_t fec = 0;
    uint64_t iters;
    size_t bec_tmp;

    FILE* fp = fopen(logfile, "w");
    if(!fp)
    {
        std::cout << "can not open logfile " << logfile << " for writing" << "\n";
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
		auto time_start = std::chrono::high_resolution_clock::now();

		do
		{
			encode_all0(x, c);

			simulate_awgn(x, y, sigma2);

			//puncturing & shortening
			if(mLdpcCode->num_puncture() != 0)
			{
				for(size_t j = 0; j < mLdpcCode->num_puncture(); j++) {
					mLdpcDecoder->mLLRIn[mLdpcCode->puncture()[j]] = 0;
				}
			}
			if(mLdpcCode->num_shorten() != 0)
			{
				for(size_t j = 0; j < mLdpcCode->num_shorten(); j++) {
					mLdpcDecoder->mLLRIn[mLdpcCode->shorten()[j]] = 99999.9;
				}
			}

			for(size_t j = 0; j < n; j++)
			{
				calc_llrs(y[j], sigma2, l_tmp);
				for(size_t k = 0; k < bits; k++) {
					mLdpcDecoder->mLLRIn[bit_mapper[k][j]] = l_tmp[k];
				}
			}

			for(size_t j = 0; j < mLdpcCode->nc(); j++) {
				mLdpcDecoder->mLLRIn[j] *= (1-2*c[j]);
			}

			//decode
			//iters += mLdpcDecoder->decode_legacy();
			//iters += mLdpcDecoder->decode_layered_legacy();
			iters += mLdpcDecoder->decode_layered();


			frames++;

			bec_tmp = 0;
			for(size_t j = 0; j < mLdpcCode->nc(); j++)
			{
				bec_tmp += (mLdpcDecoder->mLLROut[j] <= 0);
			}

			if (bec_tmp > 0)
			{
				bec += bec_tmp;
				fec++;

				auto time_dur = std::chrono::high_resolution_clock::now() - time_start;
				uint64_t t = static_cast<uint64_t>(std::chrono::duration_cast<std::chrono::microseconds>(time_dur).count());
				printf("FRAME ERROR (%lu/%lu) in frame %lu @SNR = %.3f: BER=%.2e, FER=%.2e, TIME/FRAME=%.3f ms, AVGITERS=%.2f\n",
					fec, min_fec, frames, 10*log10(1/sigma2),
					(double) bec/(frames*mLdpcCode->nc()), (double) fec/frames,
					(double) t/frames * 1e-3,
					(double) iters/frames
				);

				log_error(c, frames, snrs[i]);
			}


		} while (fec < min_fec && frames < max_frames); //end while

		fprintf(fp, "%lf %.3e %.3e %lu %.3e\n", snrs[i], (double) fec/frames, (double) bec/(frames*mLdpcCode->nc()), frames, (double) iters/frames);
		fflush(fp);
	}//end for

	fclose(fp);

	delete[] x;
	delete[] y;
	delete[] c;
	delete[] l_tmp;
}


void ldpc_sim::log_error(bits_t* c, const uint64_t frame_num, const double snr)
{
	char errors_file[MAX_FILENAME_LEN];
    snprintf(errors_file, MAX_FILENAME_LEN, "errors_%s", logfile);

    FILE *fp = fopen(errors_file, "a+");
    if(!fp)
	{
        printf("can not open error log file.\n");
        exit(EXIT_FAILURE);
    }

    /* calculation of syndrome and failed syndrome checks */
    size_t synd_weight = 0;
    for (size_t i = 0; i < mLdpcCode->mc(); i++) {
        synd_weight += (size_t) mLdpcDecoder->mSynd[i];
    }
	std::vector<size_t> failed_checks_idx(synd_weight);
    size_t j = 0;
    for (size_t i = 0; i < mLdpcCode->mc(); i++)
	{
        if(mLdpcDecoder->mSynd[i] == 1)
		{
            failed_checks_idx[j++] = i;
        }
    }

    /* calculation of failed codeword bits */
    size_t cw_dis = 0;
    for (size_t i = 0; i < mLdpcCode->nc(); i++)
	{
        #ifdef ENCODE
        cw_dis += ((mLdpcDecoder->mLLROut[i] <= 0) != c[i]);
        #else
        cw_dis += ((mLdpcDecoder->mLLROut[i] <= 0) != 0);
        #endif
    }

	size_t* x = new size_t[n];
	size_t* xhat = new size_t[n];
	bits_t* chat = new bits_t[mLdpcCode->nc()];
	for (size_t i = 0; i < mLdpcCode->nc(); i++) {
		chat[i] = (bits_t) (mLdpcDecoder->mLLROut[i] <= 0);
	}

    map_c_to_x(c, x);
    map_c_to_x(chat, xhat);
    double cw_dis_euc = 0;
    for (size_t i = 0; i < n; i++) {
        #ifdef ENCODE
        cw_dis_euc += (cstll->X[x[i]] - cstll->X[xhat[i]]) * (cstll->X[x[i]] - cstll->X[xhat[i]]);
        #else
        cw_dis_euc += (cstll->X[0] - cstll->X[xhat[i]]) * (cstll->X[0] - cstll->X[xhat[i]]);
        #endif
    }
	std::vector<size_t> failed_bits_idx(cw_dis);
    j = 0;
    for (size_t i = 0; i < mLdpcCode->nc(); i++) {
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
	for (auto failed_bits_idx_i : failed_bits_idx) {
		fprintf(fp, "%lu ", failed_bits_idx_i);
	}
    fprintf(fp, " -- ");
    fprintf(fp, "synd weight: %lu | ", synd_weight);
	for (auto failed_checks_idx_i : failed_checks_idx) {
		fprintf(fp, "%lu ", failed_checks_idx_i);
	}
    fprintf(fp, "\n");
    fclose(fp);

	delete[] x;
	delete[] xhat;
	delete[] chat;
}
