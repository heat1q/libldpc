#include "ldpcsim.h"

namespace ldpc
{
//Init constructor
__host__ constellation::constellation(labels_t pM)
    : mM(pM), mLog2M(log2(pM)), mPX(pM, 1.0 / pM), mX(pM, 0.0)
{
    double m = 0;
    for (std::size_t j = 0; j < mM; ++j)
    {
        mX[j] = (double)-mM + 1 + 2 * j;
        m += mX[j] * mX[j] * mPX[j];
    }

    for (std::size_t j = 0; j < mM; ++j)
    {
        mX[j] = mX[j] / sqrt(m);
    }
}

/*
 * ldpc_sim_device
 */
//init constructor
__host__ ldpc_sim_device::ldpc_sim_device(cuda_ptr<ldpc_code_device> &pCode, const char *pSimFileName, const char *pMapFileName, labels_t pNumThreads)
    : mLdpcCode(pCode), mLdpcDecoderVec(), mThreads(pNumThreads)
{
    try
    {
        FILE *fpSim;
        FILE *fpMap;

        std::ifstream fsSim(pSimFileName);
        std::ifstream fsMap(pMapFileName);
        std::string fsLine;
        std::string fsSubStr;

        char tmp[1000];
        double m;
        char *tokp;

        std::getline(fsSim, fsLine);
        mLogfile = fsLine.substr(fsLine.find(":") + 2); //name

        labels_t M;
        std::getline(fsSim, fsLine);
        fsSubStr = fsLine.substr(fsLine.find(":") + 2); //M
        M = static_cast<labels_t>(std::stoul(fsSubStr));

        //setup constellation
        mConstellation = constellation(M);

        std::getline(fsSim, fsLine);
        fsSubStr = fsLine.substr(fsLine.find(":") + 2); //bits
        mBits = static_cast<labels_t>(std::stoul(fsSubStr));

        mLabels = cuda_vector<labels_t>();
        std::getline(fsSim, fsLine);
        fsSubStr = fsLine.substr(fsLine.find(":") + 2); //labels
        strcpy(tmp, fsSubStr.c_str());
        tokp = strtok(tmp, ", ");
        std::size_t i = 0;
        while (tokp != nullptr)
        {
            mLabels.push_back(static_cast<labels_t>(atoi(tokp)));
            tokp = strtok(nullptr, ", ");
            i++;
        }
        if (i != M)
        {
            throw std::runtime_error("error parsing simfile: number of constellation points does not match label size");
        }

        //reverse labels
        mLabelsRev = cuda_vector<labels_t>(M, 0);
        for (std::size_t i = 0; i < M; ++i)
        {
            mLabelsRev[mLabels[i]] = i;
        }

        mSnrs = cuda_vector<double>();
        std::getline(fsSim, fsLine);
        fsSubStr = fsLine.substr(fsLine.find(":") + 2); //snrs
        strcpy(tmp, fsSubStr.c_str());
        tokp = strtok(tmp, ",");
        i = 0;
        while (tokp != nullptr)
        {
            mSnrs.push_back(static_cast<double>(atof(tokp)));
            tokp = strtok(nullptr, ",");
            i++;
        }

        std::getline(fsSim, fsLine);
        fsSubStr = fsLine.substr(fsLine.find(":") + 2); //max frames
        mMaxFrames = std::stoul(fsSubStr);

        std::getline(fsSim, fsLine);
        fsSubStr = fsLine.substr(fsLine.find(":") + 2); //min fec
        mMinFec = std::stoul(fsSubStr);

        std::getline(fsSim, fsLine);
        fsSubStr = fsLine.substr(fsLine.find(":") + 2); //bp iter
        mBPIter = std::stoul(fsSubStr);

        bool earlyTerm;
        std::getline(fsSim, fsLine);
        fsSubStr = fsLine.substr(fsLine.find(":") + 2); //early term
        earlyTerm = static_cast<bits_t>(std::stoul(fsSubStr));

        if (mLdpcCode->nct() % mBits != 0)
        {
            throw std::runtime_error("Chosen setting m with n_c does not work. Please correct.");
        }

        mN = mLdpcCode->nct() / mBits;
        mSE = (((double)mLdpcCode->kct()) / mLdpcCode->nct()) * mBits;

        //setup bitmapper
        mBitMapper = cuda_vector<cuda_vector<std::size_t>>(mBits, cuda_vector<std::size_t>(mN, 0));
        std::getline(fsMap, fsLine);
        std::size_t pos = 0;
        for (auto &ibm : mBitMapper)
        {
            for (auto &jbm : ibm)
            {
                pos = fsLine.find(",");
                fsSubStr = fsLine.substr(0, pos);
                jbm = std::stoul(fsSubStr);
                fsLine.erase(0, pos + 2); //+2 for comma & space
            }
        }

        mBitPos = cuda_vector<std::size_t>(mLdpcCode->nct(), 0);
        bits_t found_p = 0;
        bits_t found_s = 0;

        std::size_t idx = 0;
        for (std::size_t i = 0; i < mLdpcCode->nc(); i++)
        {
            for (std::size_t j = 0; j < mLdpcCode->num_shorten(); j++)
            {
                if (mLdpcCode->shorten()[j] == i)
                {
                    found_s = 1;
                    break;
                }
            }
            for (std::size_t j = 0; j < mLdpcCode->num_puncture(); j++)
            {
                if (mLdpcCode->puncture()[j] == i)
                {
                    found_p = 1;
                    break;
                }
            }
            if ((!found_s) && (!found_p))
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
        for (std::size_t i = 0; i < mThreads; i++)
        {
            mLdpcDecoderVec.push_back(cuda_ptr<ldpc_decoder_device>(ldpc_decoder_device(mLdpcCode, mBPIter, earlyTerm)));
        }

        //channel i/o
        mX = mat_size_t(mThreads, vec_size_t(mN));
        mY = mat_double_t(mThreads, vec_double_t(mN));
        mC = mat_bits_t(mThreads, vec_bits_t(mLdpcCode->nc()));

        //rng states
        mCurandStateEncoding = mat_curandState_t(mThreads, cuda_vector<curandState_t>(mLdpcCode->nct()));
        mCurandState = mat_curandState_t(mThreads, cuda_vector<curandState_t>(mN));

        //mem_prefetch();
    }
    catch (std::exception &e)
    {
        std::cout << "Error: ldpc_sim_device::ldpc_sim_device() " << e.what() << "\n";
        exit(EXIT_FAILURE);
    }
}


__host__ double ldpc_sim_device::randn()
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

__host__ __device__ double ldpc_sim_device::simulate_awgn(double pSigma2)
{
    double a = 0;
    double Pn = 0;
    double Px = 0;

    for (std::size_t i = 0; i < mN; i++)
    {
        a = randn() * sqrt(pSigma2);
        Pn += a * a;
        Px += mConstellation.X()[mX[0][i]] * mConstellation.X()[mX[0][i]];
        mY[0][i] = mConstellation.X()[mX[0][i]] + a;
    }

    return Px / Pn;
}

__host__ __device__ void ldpc_sim_device::encode_all0()
{
    for (std::size_t i = 0; i < mLdpcCode->nct(); i++)
    {
        mC[0][mBitPos[i]] = rand() & 1;
    }

    for (std::size_t i = 0; i < mLdpcCode->num_puncture(); i++)
    {
        mC[0][mLdpcCode->puncture()[i]] = rand() & 1;
    }

    for (std::size_t i = 0; i < mLdpcCode->num_shorten(); i++)
    {
        mC[0][mLdpcCode->shorten()[i]] = 0;
    }

    map_c_to_x();
}

__host__ __device__ void ldpc_sim_device::map_c_to_x()
{
    std::size_t tmp;

    for (std::size_t i = 0; i < mN; i++)
    {
        tmp = 0;
        for (std::size_t j = 0; j < mBits; j++)
        {
            tmp += mC[0][mBitMapper[j][i]] << (mBits - 1 - j);
        }

        mX[0][i] = mLabelsRev[tmp];
    }
}

__host__ __device__ void ldpc_sim_device::calc_llrs(double sigma2)
{
    std::vector<double> llr_tmp(mBits);

    for (std::size_t l = 0; l < mN; l++)
    {
        double tmp0, tmp1;

        for (std::size_t i = 0; i < mConstellation.log2M(); i++)
        {
            tmp0 = 0.0;
            tmp1 = 0.0;
            for (std::size_t j = 0; j < mConstellation.M(); j++)
            {
                if (mLabels[j] & (1 << (mConstellation.log2M() - 1 - i)))
                {
                    tmp1 += exp(-(mY[0][l] - mConstellation.X()[j]) * (mY[0][l] - mConstellation.X()[j]) / (2 * sigma2)) * mConstellation.pX()[j];
                }
                else
                {
                    tmp0 += exp(-(mY[0][l] - mConstellation.X()[j]) * (mY[0][l] - mConstellation.X()[j]) / (2 * sigma2)) * mConstellation.pX()[j];
                }
            }
            double val = log(tmp0 / tmp1);
            // check usually required when PAS is used with large constellations
            // and severely shaped distributions
            if (std::isinf(val) == +1)
            {
                llr_tmp[i] = MAX_LLR;
            }
            else if (std::isinf(val) == -1)
            {
                llr_tmp[i] = MIN_LLR;
            }
            else
            {
                llr_tmp[i] = val;
            }
        }

        for (std::size_t k = 0; k < mBits; k++)
        {
            mLdpcDecoderVec[0]->mLLRIn[mBitMapper[k][l]] = llr_tmp[k];
        }
    }

    for (std::size_t j = 0; j < mLdpcCode->nc(); j++)
    {
        mLdpcDecoderVec[0]->mLLRIn[j] *= (1 - 2 * mC[0][j]);
    }
}

__host__ void ldpc_sim_device::print()
{
    std::cout << "logfile: " << mLogfile << "\n";
    printf("n: %lu\n", mN);
    printf("M: %hu\n", mConstellation.M());
    printf("bits: %hu\n", mBits);
    printf("SNRs: ");
    for (const auto &i : mSnrs)
    {
        printf("%lf ", i);
    }
    printf("\n");
    printf("labels:\n");
    for (std::size_t i = 0; i < mConstellation.M(); i++)
    {
        printf("\t%lu: ", i);
        dec2bin(mLabels[i], mBits);
        printf("\n");
    }
    printf("max frames: %lu\n", mMaxFrames);
    printf("min fec: %lu\n", mMinFec);
    printf("bp iter: %lu\n", mBPIter);
    printf("early term: %u\n", mLdpcDecoderVec[0]->early_termination());
    printf("SE: %.4lf\n", mSE);
}

__host__ void ldpc_sim_device::log_error(std::size_t pFrameNum, double pSNR, labels_t pThreads)
{
    char errors_file[MAX_FILENAME_LEN];
    snprintf(errors_file, MAX_FILENAME_LEN, "errors_%s", mLogfile.c_str());

    FILE *fp = fopen(errors_file, "a+");
    if (!fp)
    {
        printf("can not open error log file.\n");
        exit(EXIT_FAILURE);
    }

    // calculation of syndrome and failed syndrome checks
    std::size_t synd_weight = 0;
    for (auto si : mLdpcDecoderVec[pThreads]->mSynd)
    {
        synd_weight += si;
    }
    std::vector<std::size_t> failed_checks_idx(synd_weight);
    std::size_t j = 0;
    for (std::size_t i = 0; i < mLdpcCode->mc(); i++)
    {
        if (mLdpcDecoderVec[pThreads]->mSynd[i] == 1)
        {
            failed_checks_idx[j++] = i;
        }
    }

    // calculation of failed codeword bits
    std::size_t cw_dis = 0;
    for (std::size_t i = 0; i < mLdpcCode->nc(); i++)
    {
#ifdef ENCODE
        cw_dis += ((mLdpcDecoderVec[pThreads]->mLLROut[i] <= 0) != mC[pThreads][i]);
#else
        cw_dis += ((mLdpcDecoderVec[pThreads]->mLLROut[i] <= 0) != 0);
#endif
    }

    std::vector<std::size_t> x(mN);
    std::vector<std::size_t> xhat(mN);
    std::vector<std::size_t> chat(mLdpcCode->nc());

    for (std::size_t i = 0; i < mLdpcCode->nc(); ++i)
    {
        chat[i] = (mLdpcDecoderVec[pThreads]->mLLROut[i] <= 0);
    }

    std::size_t tmp;
    //map c to x map_c_to_x(c, x);
    for (std::size_t i = 0; i < mN; i++)
    {
        tmp = 0;
        for (std::size_t j = 0; j < mBits; j++)
        {
            tmp += mC[pThreads][mBitMapper[j][i]] << (mBits - 1 - j);
        }

        x[i] = mLabelsRev[tmp];
    }

    //map_c_to_x(chat, xhat);
    for (std::size_t i = 0; i < mN; i++)
    {
        tmp = 0;
        for (std::size_t j = 0; j < mBits; j++)
        {
            tmp += chat[mBitMapper[j][i]] << (mBits - 1 - j);
        }

        xhat[i] = mLabelsRev[tmp];
    }

    double cw_dis_euc = 0;
    for (std::size_t i = 0; i < mN; i++)
    {
#ifdef ENCODE
        cw_dis_euc += (mConstellation.X()[x[i]] - mConstellation.X()[xhat[i]]) * (mConstellation.X()[x[i]] - mConstellation.X()[xhat[i]]);
#else
        cw_dis_euc += (mConstellation.X()[0] - mConstellation.X()[xhat[i]]) * (mConstellation.X()[0] - mConstellation.X()[xhat[i]]);
#endif
    }
    std::vector<std::size_t> failed_bits_idx(cw_dis);
    j = 0;
    for (std::size_t i = 0; i < mLdpcCode->nc(); i++)
    {
#ifdef ENCODE
        if (chat[i] != mC[pThreads][i])
        {
            failed_bits_idx[j++] = i;
        }
#else
        if (chat[i] != 0)
        {
            failed_bits_idx[j++] = i;
        }
#endif
    }

    // print results in file
    fprintf(fp, "SNR: %.2f -- frame: %lu -- is codeword: %d -- dE(c,chat): %.3f -- dH(c,chat): %lu | ", pSNR, pFrameNum, synd_weight == 0, cw_dis_euc, cw_dis);
    for (auto failed_bits_idx_i : failed_bits_idx)
    {
        fprintf(fp, "%lu ", failed_bits_idx_i);
    }
    fprintf(fp, " -- ");
    fprintf(fp, "synd weight: %lu | ", synd_weight);
    for (auto failed_checks_idx_i : failed_checks_idx)
    {
        fprintf(fp, "%lu ", failed_checks_idx_i);
    }
    fprintf(fp, "\n");
    fclose(fp);
}

//measure constant time w/o decoding with pCount samples
//returns time in us
#ifdef LOG_TP
__host__ std::size_t ldpc_sim_device::frame_const_time(double pSigma2, std::size_t pCount)
{
#ifdef USE_LEGACY_DEC || defined USE_CPU_FRAME //for cpu frame proc.
    auto tconstStart = std::chrono::high_resolution_clock::now();

    for (std::size_t i = 0; i < pCount; ++i)
    {
        encode_all0();
        simulate_awgn(pSigma2);
        calc_llrs(pSigma2);

        std::size_t tmp = 0;
        for (std::size_t j = 0; j < mLdpcCode->nc(); ++j)
        {
            tmp += (mLdpcDecoderVec[0]->mLLROut[j] <= 0);
        }
    }
#else
    //call one time to reduce memory overhead
    cudakernel::sim::frame_time<<<mThreads, 1>>>(this, pSigma2);
    cudaDeviceSynchronize();

    auto tconstStart = std::chrono::high_resolution_clock::now();

    for (std::size_t i = 0; i < pCount; ++i)
    {
        cudakernel::sim::frame_time<<<mThreads, 1>>>(this, pSigma2);//launch w/o decoding i.e zero iterations
	cudaDeviceSynchronize();
    }

    for (std::size_t i = 0; i < mThreads * pCount; ++i)
    {
        std::size_t tmp = 0;
        for (std::size_t j = 0; j < mLdpcCode->nc(); ++j)
        {
            tmp += (mLdpcDecoderVec[0]->mLLROut[j] <= 0);
        }
    }
#endif
    auto tconstDiff = std::chrono::high_resolution_clock::now() - tconstStart;
    std::size_t tconst = static_cast<std::size_t>(std::chrono::duration_cast<std::chrono::microseconds>(tconstDiff).count());

    return tconst / pCount;
}
#endif
} // namespace ldpc
