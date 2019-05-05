#include "../ldpcsim.h"

#include <math.h>
#include <exception>

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
            mLdpcDecoderVec.push_back(cuda_ptr<ldpc_decoder_device>(ldpc_decoder_device(mLdpcCode, mBPIter, true)));
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

__host__ void ldpc_sim_device::mem_prefetch()
{
    /*
    mConstellation.mem_prefetch();

    for (auto &bmi : mBitMapper)
    {
        bmi.mem_prefetch();
    }
    mBitMapper.mem_prefetch();

    mBitPos.mem_prefetch();
    mLabelsRev.mem_prefetch();
    mLabels.mem_prefetch();
    mSnrs.mem_prefetch();

    mX[0].mem_prefetch();
    mY[0].mem_prefetch();
    mC[0].mem_prefetch();

    mCurandStateEncoding.mem_prefetch();
    mCurandState.mem_prefetch();
*/
}

//start simulation for parallel frame processing on gpu
//specified with mThreads
__host__ void ldpc_sim_device::start_device()
{
    //setup random number generators on device
    cudakernel::sim::setup_rng<<<mThreads, NUMK_THREADS>>>(this);
    cudaDeviceSynchronize();

    double sigma2;
    std::size_t frames;
    std::size_t bec = 0;
    std::size_t fec = 0;
    std::size_t iters;
    std::size_t bec_tmp;

    std::vector<std::string> printResStr(mSnrs.size() + 1, std::string());
    std::ofstream fp;
    char resStr[128];

#ifdef LOG_FRAME_TIME
    printResStr[0].assign("snr fer ber frames avg_iter frame_time");
#elif defined LOG_TP
    printResStr[0].assign("snr fer ber frames avg_iter frame_time dec_time throughput");
#else
    printResStr[0].assign("snr fer ber frames avg_iter");
#endif

    for (std::size_t i = 0; i < mSnrs.size(); ++i)
    {
        bec = 0;
        fec = 0;
        frames = 0;
        iters = 0;
        sigma2 = pow(10, -mSnrs[i] / 10);
        auto timeStart = std::chrono::high_resolution_clock::now();

#ifdef LOG_TP
        std::size_t tconst = frame_const_time(sigma2, 10);
#endif
        do
        {
            //launch the frame processing kernel
            cudakernel::sim::frame_proc<<<mThreads, 1>>>(this, sigma2);
            cudaDeviceSynchronize();

            //now check the processed frames
            for (labels_t k = 0; k < mThreads; ++k) //TODO adjust for parallel threads
            {
                if (fec >= mMinFec)
                {
                    break;
                } //prevent checking more frames than necessary

                iters += mLdpcDecoderVec[k]->mIter;
                ++frames;

                bec_tmp = 0;
                for (std::size_t j = 0; j < mLdpcCode->nc(); ++j)
                {
                    bec_tmp += (mLdpcDecoderVec[k]->mLLROut[j] <= 0);
                }

                if (bec_tmp > 0)
                {
                    bec += bec_tmp;
                    ++fec;

                    auto timeNow = std::chrono::high_resolution_clock::now();
                    auto timeFrame = timeNow - timeStart; //eliminate const time for printing etc
                    std::size_t tFrame = static_cast<std::size_t>(std::chrono::duration_cast<std::chrono::microseconds>(timeFrame).count());
                    tFrame = tFrame / (frames + mThreads - k - 1);
#ifdef LOG_TP
                    printf("\r %2lu/%2lu  |  %12lu  |  %.3f  |  %.2e  |  %.2e  |  %.1e  |  %.3fms  |  %6us  |  %.2fMbits/s",
                           fec, mMinFec, frames, mSnrs[i],
                           (double)bec / (frames * mLdpcCode->nc()), (double)fec / frames,
                           (double)iters / frames,
                           (double)tFrame * 1e-3,
                           (double)tFrame - tconst, //decoding time
                           (double)mLdpcCode->nc() / tFrame); //decoding throughput
#else
                    printf("\r %2lu/%2lu  |  %12lu  |  %.3f  |  %.2e  |  %.2e  |  %.1e  |  %.3fms",
                           fec, mMinFec, frames, mSnrs[i],
                           (double)bec / (frames * mLdpcCode->nc()), (double)fec / frames,
                           (double)iters / frames,
                           (double)tFrame * 1e-3);
#endif
                    fflush(stdout);

#ifdef LOG_FRAME_TIME
                    sprintf(resStr, "%lf %.3e %.3e %lu %.3e %.6f",
                            mSnrs[i], (double)fec / frames, (double)bec / (frames * mLdpcCode->nc()),
                            frames, (double)iters / frames, (double)tFrame * 1e-6);
#elif defined LOG_TP
                    sprintf(resStr, "%lf %.3e %.3e %lu %.3e %.6f %.6f %lu",
                            mSnrs[i], (double)fec / frames, (double)bec / (frames * mLdpcCode->nc()),
                            frames, (double)iters / frames, (double)tFrame * 1e-6,
                            (double)tFrame - tconst * 1e-6, (double)mLdpcCode->nc() / (tFrame * 1e-6));
#else
                    sprintf(resStr, "%lf %.3e %.3e %lu %.3e",
                            mSnrs[i], (double)fec / frames, (double)bec / (frames * mLdpcCode->nc()),
                            frames, (double)iters / frames);
#endif
                    printResStr[i + 1].assign(resStr);

                    try
                    {
                        fp.open(mLogfile);
                        for (const auto &x : printResStr)
                        {
                            fp << x << "\n";
                        }
                        fp.close();
                    }
                    catch (...)
                    {
                        std::cout << "Warning: can not open logfile " << mLogfile << " for writing"
                                  << "\n";
                    }

                    log_error(frames, mSnrs[i], k);

                    timeStart += std::chrono::high_resolution_clock::now() - timeNow; //dont measure time for printing files
                }
            }
        } while (fec < mMinFec && frames < mMaxFrames); //end while
        printf("\n");
    } //end for
}

//start simulation on cpu
__host__ void ldpc_sim_device::start()
{
    double sigma2;
    std::size_t frames;
    std::size_t bec = 0;
    std::size_t fec = 0;
    std::size_t iters;
    std::size_t bec_tmp;

    std::vector<std::string> printResStr(mSnrs.size() + 1, std::string());
    std::ofstream fp;
    char resStr[128];

#ifdef LOG_FRAME_TIME
    printResStr[0].assign("snr fer ber frames avg_iter time_frame[ms]");
#else
    printResStr[0].assign("snr fer ber frames avg_iter");
#endif

    for (std::size_t i = 0; i < mSnrs.size(); ++i)
    {
        bec = 0;
        fec = 0;
        frames = 0;
        iters = 0;
        sigma2 = pow(10, -mSnrs[i] / 10);
        auto timeStart = std::chrono::high_resolution_clock::now();

#ifdef LOG_TP
        std::size_t tconst = frame_const_time(sigma2, 10);
#endif
        do
        {
            encode_all0();
            simulate_awgn(sigma2);

            //puncturing & shortening
            if (mLdpcCode->num_puncture() != 0)
            {
                for (std::size_t j = 0; j < mLdpcCode->num_puncture(); j++)
                {
                    mLdpcDecoderVec[0]->mLLRIn[mLdpcCode->puncture()[j]] = 0;
                }
            }
            if (mLdpcCode->num_shorten() != 0)
            {
                for (std::size_t j = 0; j < mLdpcCode->num_shorten(); j++)
                {
                    mLdpcDecoderVec[0]->mLLRIn[mLdpcCode->shorten()[j]] = 99999.9;
                }
            }

            calc_llrs(sigma2);

            //decode
#ifdef USE_LEGACY_DEC
            iters += mLdpcDecoderVec[0]->decode_legacy();
#else
            iters += mLdpcDecoderVec[0]->decode_layered();
#endif

            ++frames;

            bec_tmp = 0;
            for (std::size_t j = 0; j < mLdpcCode->nc(); ++j)
            {
                bec_tmp += (mLdpcDecoderVec[0]->mLLROut[j] <= 0);
            }

            if (bec_tmp > 0)
            {
                bec += bec_tmp;
                ++fec;

                auto timeNow = std::chrono::high_resolution_clock::now();
                auto timeFrame = timeNow - timeStart; //eliminate const time for printing etc
                std::size_t tFrame = static_cast<std::size_t>(std::chrono::duration_cast<std::chrono::microseconds>(timeFrame).count());
                tFrame = tFrame / frames;
#ifdef LOG_TP
                printf("\r %2lu/%2lu  |  %12lu  |  %.3f  |  %.2e  |  %.2e  |  %.1e  |  %.3fms  |  %6us  |  %.2fMbits/s",
                       fec, mMinFec, frames, mSnrs[i],
                       (double)bec / (frames * mLdpcCode->nc()), (double)fec / frames,
                       (double)iters / frames,
                       (double)tFrame * 1e-3,
                       (double)tFrame - tconst, //decoding time
                       (double)mLdpcCode->nc() / tFrame); //decoding throughput
#else
                printf("\r %2lu/%2lu  |  %12lu  |  %.3f  |  %.2e  |  %.2e  |  %.1e  |  %.3fms",
                       fec, mMinFec, frames, mSnrs[i],
                       (double)bec / (frames * mLdpcCode->nc()), (double)fec / frames,
                       (double)iters / frames,
                       (double)tFrame * 1e-3);
#endif
                fflush(stdout);

#ifdef LOG_FRAME_TIME
                sprintf(resStr, "%lf %.3e %.3e %lu %.3e %.6f",
                        mSnrs[i], (double)fec / frames, (double)bec / (frames * mLdpcCode->nc()),
                        frames, (double)iters / frames, (double)tFrame * 1e-6);
#elif defined LOG_TP
                sprintf(resStr, "%lf %.3e %.3e %lu %.3e %.6f %.6f %lu",
                        mSnrs[i], (double)fec / frames, (double)bec / (frames * mLdpcCode->nc()),
                        frames, (double)iters / frames, (double)tFrame * 1e-6,
                        (double)tFrame - tconst * 1e-6, (double)mLdpcCode->nc() / (tFrame * 1e-6));
#else
                sprintf(resStr, "%lf %.3e %.3e %lu %.3e",
                        mSnrs[i], (double)fec / frames, (double)bec / (frames * mLdpcCode->nc()),
                        frames, (double)iters / frames);
#endif
                printResStr[i + 1].assign(resStr);

                try
                {
                    fp.open(mLogfile);
                    for (const auto &x : printResStr)
                    {
                        fp << x << "\n";
                    }
                    fp.close();
                }
                catch (...)
                {
                    std::cout << "Warning: can not open logfile " << mLogfile << " for writing"
                              << "\n";
                }

                log_error(frames, mSnrs[i], 0);

                timeStart += std::chrono::high_resolution_clock::now() - timeNow; //dont measure time for printing files
            }
        } while (fec < mMinFec && frames < mMaxFrames); //end while
        printf("\n");
    } //end for
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
    printf("SE: %.4lf\n", mSE);
    //for (std::size_t i = 0; i < mBits; i++)
    //{
    //    printf("Bit Mapping B%lu: ", i);
    //    for (std::size_t j = 0; j < mN; j++)
    //    {
    //        printf("%lu, ", mBitMapper[i][j]);
    //    }
    //    printf("\n");
    //}
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

        std::size_t tmo = 0;
        for (std::size_t j = 0; j < mLdpcCode->nc(); ++j)
        {
            bec_tmp += (mLdpcDecoderVec[0]->mLLROut[j] <= 0);
        }
    }

    auto tconstDiff = std::chrono::high_resolution_clock::now() - tconstStart;
    std::size_t tconst = static_cast<std::size_t>(std::chrono::duration_cast<std::chrono::microseconds>(tconstDiff).count());

    return tconst / pCount;
#else
    //call one time to reduce memory overhead
    cudakernel::sim::frame_time<<<mThreads, 1>>>(this, pSigma2, 1);
    cudaDeviceSynchronize();

    auto tconstStart = std::chrono::high_resolution_clock::now();

    cudakernel::sim::frame_time<<<mThreads, 1>>>(this, pSigma2, pCount);
    cudaDeviceSynchronize();

    auto tconstDiff = std::chrono::high_resolution_clock::now() - tconstStart;
    std::size_t tconst = static_cast<std::size_t>(std::chrono::duration_cast<std::chrono::microseconds>(tconstDiff).count());

    return tconst / pCount / mThreads;
#endif
}
#endif
} // namespace ldpc
