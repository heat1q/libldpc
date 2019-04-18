#include "../ldpcsim.h"

#include <math.h>
#include <exception>

using namespace ldpc;

//Init constructor
__host__ constellation::constellation(const uint16_t pM)
	: mM(pM), mLog2M(log2(pM))
	, mPX(pM, 1.0/pM), mX(pM, 0.0)
{
	double m = 0;
	for (size_t j = 0; j < mM; ++j)
	{
		mX[j] = (double) -mM+1+2*j;
		m += mX[j] * mX[j] * mPX[j];
	}

	for (size_t j = 0; j < mM; ++j) {
		mX[j] = mX[j]/sqrt(m);
	}
}


/*
 * ldpc_sim_device
*/
//init constructor
ldpc_sim_device::ldpc_sim_device(cudamgd_ptr<ldpc_code_device> pCode, const char* pSimFileName, const char* pMapFileName)
: mLdpcCode(pCode), mLdpcDecoder(nullptr)
{
	try
	{
		FILE *fpSim;
		FILE *fpMap;

		char tmp[1000];
		double m;
		char* tokp;

		fpSim = fopen(pSimFileName, "r");
		if(!fpSim) { throw std::runtime_error("can not open sim file"); }

		fpMap = fopen(pMapFileName, "r");
		if(!fpMap) { throw std::runtime_error("can not open mapping file"); }

		fscanf(fpSim, "name: %256s\n", mLogfile);

		uint16_t M;
		fscanf(fpSim, "M: %hu\n", &M);
		mConstellation = constellation(M); //setup constellation
		fscanf(fpSim, "bits: %hu\n", &mBits);

		mLabels = vector_mgd<uint16_t>();
		fscanf(fpSim, "labels: %[^\n]\n", tmp);
		tokp = strtok(tmp, ", ");
		size_t i = 0;
		while(tokp != nullptr)
		{
			mLabels.push_back((unsigned short) atoi(tokp));
			tokp = strtok(nullptr, ", ");
			i++;
		}
		if(i != M)
		{
			throw std::runtime_error("error parsing simfile: number of constellation points does not match label size");
		}

		//reverse labels
		mLabelsRev = vector_mgd<uint16_t>(M, 0);
		for (size_t i = 0; i < M; ++i) {
			mLabelsRev[mLabels[i]] = i;
		}

		mSnrs = vector_mgd<double>();
		fscanf(fpSim, "snrs: %[^\n]\n", tmp);
		tokp = strtok(tmp, ",");
		i = 0;
		while(tokp != nullptr)
		{
			mSnrs.push_back(static_cast<double> (atof(tokp)));
			tokp = strtok(nullptr, ",");
			i++;
		}

		fscanf(fpSim, "max frames: %lu\n", &mMaxFrames);
		fscanf(fpSim, "min fec: %lu\n", &mMinFec);
		fscanf(fpSim, "bp iter: %lu\n", &mBPIter);

		if(mLdpcCode->nct() % mBits != 0)
		{
			throw std::runtime_error("Chosen setting m with n_c does not work. Please correct.");
		}

		mN = mLdpcCode->nct()/mBits;
		mSE = (((double) mLdpcCode->kct())/mLdpcCode->nct()) * mBits;


		mBitMapper = vector_mgd< vector_mgd<size_t> >(mBits, vector_mgd<size_t>(mN, 0));
		for (auto& ibm : mBitMapper)
		{
			for (auto& jbm : ibm) {
				fscanf(fpMap, " %lu, ", &jbm);
			}
		}


		mBitPos = vector_mgd<size_t>(mLdpcCode->nct(), 0);
		bits_t found_p = 0;
		bits_t found_s = 0;

		size_t idx = 0;
		for(size_t i = 0; i < mLdpcCode->nc(); i++)
		{
			for(size_t j = 0; j < mLdpcCode->num_shorten(); j++)
			{
				if(mLdpcCode->shorten()[j] == i)
				{
					found_s = 1;
					break;
				}
			}
			for(size_t j = 0; j < mLdpcCode->num_puncture(); j++)
			{
				if(mLdpcCode->puncture()[j] == i)
				{
					found_p = 1;
					break;
				}
			}
			if((!found_s) && (!found_p))
			{
				mBitPos[idx] = i;
				idx++;
			}
			else
			{
				found_p = 0;
				found_s = 0;
			}
		}

		//changed with each frame
		//set up decoder
		mLdpcDecoder = cudamgd_ptr<ldpc_decoder_device>(
			ldpc_decoder_device(mLdpcCode, mBPIter, true)
		);

		//channel i/o
		mX = vec_size_t(mN);
		mY = vec_double_t(mN);
		mC = vec_bits_t(mLdpcCode->nc());
		mLTmp = vec_double_t(mBits);

		fclose(fpSim);
		fclose(fpMap);
	}
	catch(std::exception &e)
	{
		std::cout << "Error: " << e.what() << "\n";
		exit(EXIT_FAILURE);
	}
}


void ldpc_sim_device::start()
{
    double sigma2;
    uint64_t frames;
    uint64_t bec = 0;
    uint64_t fec = 0;
    uint64_t iters;
    size_t bec_tmp;

	std::vector<std::string> printResStr(mSnrs.size()+1, std::string());
	std::ofstream fp;
	char resStr[128];

	#ifdef LOG_FRAME_TIME
	printResStr[0].assign("snr fer ber frames avg_iter time_frame[ms]");
	#else
	printResStr[0].assign("snr fer ber frames avg_iter");
	#endif

	/*
	* START: SIMULATION PART
	*/

	for (size_t i = 0; i < mSnrs.size(); ++i)
	{
		bec = 0;
		fec = 0;
		frames = 0;
		iters = 0;
		sigma2 = pow(10, -mSnrs[i]/10);
		auto time_start = std::chrono::high_resolution_clock::now();

		do
		{
			encode_all0(); //40ms
			simulate_awgn(sigma2); //100ms

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

			calc_llrs(sigma2); //150ms

			//decode
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
				size_t t = static_cast<size_t>(std::chrono::duration_cast<std::chrono::microseconds>(time_dur).count());
				printf("FRAME ERROR (%lu/%lu) in frame %lu @SNR = %.3f: BER=%.2e, FER=%.2e, TIME/FRAME=%.3fms, AVGITERS=%.2f\n",
					fec, mMinFec, frames, mSnrs[i],
					(double) bec/(frames*mLdpcCode->nc()), (double) fec/frames,
					(double) t/frames * 1e-3,
					(double) iters/frames
				);

				#ifdef LOG_FRAME_TIME
				sprintf(resStr, "%lf %.3e %.3e %lu %.3e %.3f"
					, mSnrs[i]
					, (double) fec/frames
					, (double) bec/(frames*mLdpcCode->nc())
					, frames, (double) iters/frames
					, (double) t/frames * 1e-3
				);
				#else
				sprintf(resStr, "%lf %.3e %.3e %lu %.3e"
					, mSnrs[i]
					, (double) fec/frames
					, (double) bec/(frames*mLdpcCode->nc())
					, frames, (double) iters/frames
				);
				#endif

				printResStr[i+1].assign(resStr);

				try
				{
					fp.open(mLogfile);
					for (const auto& x : printResStr)
					{
						fp << x << "\n";
					}
					fp.close();
				}
				catch(...)
				{
					std::cout << "Warning: can not open logfile " << mLogfile << " for writing" << "\n";
				}

				//log_error(c, frames, snrs[i]);
			}
		} while (fec < mMinFec && frames < mMaxFrames); //end while
	}//end for
}


double ldpc_sim_device::randn()
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


double ldpc_sim_device::simulate_awgn(double sigma2)
{
    double a = 0;
    double Pn = 0;
    double Px = 0;

    for (size_t i = 0; i < mN; i++)
    {
        a = randn() * sqrt(sigma2);
        Pn += a * a;
        Px += mConstellation.X()[mX[i]] * mConstellation.X()[mX[i]];
        mY[i] = mConstellation.X()[mX[i]] + a;
    }

    return Px/Pn;
}


void ldpc_sim_device::encode_all0()
{
    for(size_t i = 0; i < mLdpcCode->nct(); i++) {
        mC[mBitPos[i]] = rand() & 1;
	}

    for(size_t i = 0; i < mLdpcCode->num_puncture(); i++) {
        mC[mLdpcCode->puncture()[i]] = rand() & 1;
	}

    for(size_t i = 0; i < mLdpcCode->num_shorten(); i++) {
        mC[mLdpcCode->shorten()[i]] = 0;
	}

    map_c_to_x();
}


void ldpc_sim_device::map_c_to_x()
{
    size_t tmp;

    for(size_t i = 0; i < mN; i++)
    {
        tmp = 0;
        for(size_t j = 0; j < mBits; j++) {
            tmp += mC[mBitMapper[j][i]] << (mBits-1-j);
		}

        mX[i] = mLabelsRev[tmp];
    }
}


void ldpc_sim_device::calc_llrs(double sigma2)
{
	for(size_t l = 0; l < mN; l++)
	{
		double tmp0, tmp1;

		for(size_t i = 0; i < mConstellation.log2M(); i++)
		{
			tmp0 = 0.0;
			tmp1 = 0.0;
			for(size_t j = 0; j < mConstellation.M(); j++) {
				if(mLabels[j] & (1 << (mConstellation.log2M()-1-i))) {
					tmp1 += exp(-(mY[l]-mConstellation.X()[j])*(mY[l]-mConstellation.X()[j])/(2*sigma2)) * mConstellation.pX()[j];
				} else {
					tmp0 += exp(-(mY[l]-mConstellation.X()[j])*(mY[l]-mConstellation.X()[j])/(2*sigma2)) * mConstellation.pX()[j];
				}
			}
			double val = log(tmp0/tmp1);
			// check usually required when PAS is used with large constellations
			// and severely shaped distributions
			if(std::isinf(val) == +1) {
				mLTmp[i] = MAX_LLR;
			} else if(std::isinf(val) == -1) {
				mLTmp[i] = MIN_LLR;
			} else {
				mLTmp[i] = val;
			}
		}


		for(size_t k = 0; k < mBits; k++) {
			mLdpcDecoder->mLLRIn[mBitMapper[k][l]] = mLTmp[k];
		}
	}

	for(size_t j = 0; j < mLdpcCode->nc(); j++) {
		mLdpcDecoder->mLLRIn[j] *= (1-2*mC[j]);
	}
}


void ldpc_sim_device::print()
{
    mLdpcCode->print();

    printf("=========== SIM ===========\n");
    printf("logfile: %s\n", mLogfile);
    printf("n: %lu\n", mN);
    printf("M: %hu\n", mConstellation.M());
    printf("bits: %hu\n", mBits);
    printf("SNRs: ");
	for (const auto& i : mSnrs) {
		printf("%lf ", i);
	}
    printf("\n");
    printf("labels:\n");
    for(size_t i = 0; i < mConstellation.M(); i++)
    {
        printf("\t%lu: ", i);
        dec2bin(mLabels[i], mBits);
        printf("\n");
    }
    printf("max frames: %lu\n", mMaxFrames);
    printf("min fec: %lu\n", mMinFec);
    printf("bp iter: %lu\n", mBPIter);
    printf("SE: %.4lf\n", mSE);
    // for(size_t i = 0; i < sim.bits; i++) {
    //     printf("Bit Mapping B%lu: ", i);
    //     for(size_t j = 0; j < sim.n; j++) {
    //         printf("%lu, ", sim.bit_mapper[i][j]);
    //     }
    //     printf("\n");
    // }
    printf("=========== SIM: END ===========\n");
}
