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


//copy constructor
__host__ __device__ constellation::constellation(const constellation& pCopy)
	: mM(pCopy.mM), mLog2M(pCopy.mLog2M)
	, mPX(pCopy.mPX), mX(pCopy.mX) {}



/*
 * ldpc_sim_device
*/
//init constructor
ldpc_sim_device::ldpc_sim_device(cudamgd_ptr<ldpc_code> pCode, const char* pSimFileName, const char* pMapFileName)
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

		//set up decoder
		//mLdpcDecoder = ldpc_decoder(mLdpcCode, mBPIter, true);

		fclose(fpSim);
		fclose(fpMap);
	}
	catch(std::exception &e)
	{
		std::cout << "Error: " << e.what() << "\n";
		exit(EXIT_FAILURE);
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
