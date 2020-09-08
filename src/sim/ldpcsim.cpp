#include "ldpcsim.h"
#include "../core/functions.h"

namespace ldpc
{


ldpc_sim::ldpc_sim(const ldpc_code *code,
                   const std::string &outFile,
                   const vec_double_t &snrRange,
                   const unsigned numThreads,
                   const u64 seed,
                   const enum channel_type channel,
                   const unsigned iters,
                   const u64 maxFrames,
                   const u64 fec,
                   const bool earlyTerm)
    : ldpc_sim(code, outFile, snrRange, numThreads, seed, channel, iters, maxFrames, fec, earlyTerm, nullptr) {}

//init constructor
ldpc_sim::ldpc_sim(const ldpc_code *code,
             const std::string &outFile,
             const vec_double_t &snrRange,
             const unsigned numThreads,
             const u64 seed,
             const enum channel_type channel,
             const unsigned iters,
             const u64 maxFrames,
             const u64 fec,
             const bool earlyTerm,
             sim_results_t *results)
    : mLdpcCode(code), mLogfile(outFile), mThreads(numThreads), 
      mBPIter(iters), mMaxFrames(maxFrames), mMinFec(fec), mChannel(channel),
      mX(numThreads, vec_double_t(code->nct(), 1.0)), mY(numThreads, vec_double_t(code->nct())), //for many threads we need independent vectors
      mC(numThreads, vec_bits_t(code->nc())), mRNGSeed(seed), mResults(results)
{
    try
    {
        // initialize snr vector
        double val = snrRange[0];
        while (val < snrRange[1])
        {
            mSnrs.push_back(val);
            val += snrRange[2];
        }

        // position of transmitted bits
        for (u64 i = 0; i < mLdpcCode->nc(); i++)
        {
            auto tmp = std::find(mLdpcCode->shorten().cbegin(), mLdpcCode->shorten().cend(), i);
            if (tmp != mLdpcCode->shorten().cend()) continue; // skip if current index shortened
            tmp = std::find(mLdpcCode->puncture().cbegin(), mLdpcCode->puncture().cend(), i);
            if (tmp != mLdpcCode->puncture().cend()) continue; // skip if current index punctured

            mBitPos.push_back(i);
        }

        // RNG setup
        //std::random_device rd;
        //results may vary with same seed, since some threads are executed more than others
        for (unsigned i = 0; i < mThreads; ++i)
        {
            //decoder
            mLdpcDecoder.push_back(ldpc_decoder(mLdpcCode, this, mBPIter, earlyTerm));
            mRNG.push_back(std::mt19937(mRNGSeed + i)); // different seeds for threads
            mRandNormal.push_back(std::normal_distribution<double>(0.0, 1.0));
        }
    }
    catch (std::exception &e)
    {
        std::cout << "Error: ldpc_sim::ldpc_sim() " << e.what() << "\n";
        exit(EXIT_FAILURE);
    }
}

void ldpc_sim::simulate_awgn(double pSigma2, unsigned threadid)
{
    //double a = 0;
    //double Pn = 0;
    //double Px = 0;

    for (u64 i = 0; i < mLdpcCode->nct(); i++)
    {
        //Pn += a * a;
        //Px += mConstellation.X()[mX[i]] * mConstellation.X()[mX[i]];
        mY[threadid][i] = mRandNormal[threadid](mRNG[threadid]) * sqrt(pSigma2) + mX[threadid][i];
    }
}


std::ostream &operator<<(std::ostream &os, const ldpc_sim &sim)
{
    os << "result output file: " << sim.mLogfile << "\n";
    os << "threads: " << sim.mThreads << "\n";
    os << "snrs: ";
    for (auto x: sim.mSnrs) os << x << ", ";
    os << "\n";
    os << "max frames: " << sim.mMaxFrames << "\n";
    os << "min fec: " << sim.mMinFec << "\n";
    os << "iterations: " << sim.mBPIter <<"\n";
    os << "RNG: mt19937\n";
    os << " Thread ID | Seed\n";
    for (unsigned i = 0; i < sim.mThreads; i++)
    {
        os << " " << i << "         | " << sim. mRNGSeed + i<< "\n";
    }
    return os;
}


/*
void ldpc_sim::print_file_header(const char *binaryFile, const char *codeFile, const char *simFile, const char *mapFile)
{
    FILE *fp = fopen(mLogfile, "a+");
    fprintf(fp, "%% binary: %s (Version: %s, Built: %s)\n", binaryFile, VERSION, BUILD_DATE);
    fprintf(fp, "%% sim file: %s\n", simFile);
    fprintf(fp, "%% code file: %s\n", codeFile);
    fprintf(fp, "%% mapping file: %s\n", mapFile);
    fprintf(fp, "%% result file: %s\n", mLogfile);
    fprintf(fp, "%% iter: %lu\n", mBPIter);
    fprintf(fp, "%% max frames: %lu\n", mMaxFrames);
    fprintf(fp, "%% min fec: %lu\n", mMinFec);
    fprintf(fp, "%% BP early terminate: %hu\n", 1);
    fprintf(fp, "%% num threads: %d\n", 1);
}
*/
/*
void ldpc_sim::log_error(u64 pFrameNum, double pSNR)
{

    char errors_file[MAX_FILENAME_LEN];
    snprintf(errors_file, MAX_FILENAME_LEN, "errors_%s", mLogfile.c_str());

    FILE *fp = fopen(errors_file, "a+");
    if (!fp)


    {


        printf("can not

 open error log file.\n");
        exit(EXIT_FAILU

RE);
    }

    // calculation of syndrome and failed syndrome checks
    u64 synd_weight = 0;
    for (auto si : mLdpcDecoder->syndrome())
    {
        synd_weight += si;
    }
    std::vector<u64> failed_checks_idx(synd_weight);
    u64 j = 0;
    for (u64 i = 0; i < mLdpcCode->mc(); i++)
    {
        if (mLdpcDecoder->syndrome()[i] == 1)
        {
            failed_checks_idx[j++] = i;
        }
    }

    // calculation of failed codeword bits
    u64 cw_dis = 0;
    for (u64 i = 0; i < mLdpcCode->nc(); i++)
    {
#ifdef ENCODE
        cw_dis += ((mLdpcDecoder->llr_out()[i] <= 0) != mC[pThreads][i]);
#else
        cw_dis += ((mLdpcDecoder->llr_out()[i] <= 0) != 0);
#endif
    }

    std::vector<u64> x(mN);
    std::vector<u64> xhat(mN);
    std::vector<u64> chat(mLdpcCode->nc());

    for (u64 i = 0; i < mLdpcCode->nc(); ++i)
    {
        chat[i] = (mLdpcDecoder->llr_out()[i] <= 0);
    }

    u64 tmp;
    //map c to x map_c_to_x(c, x);
    for (u64 i = 0; i < mN; i++)
    {
        tmp = 0;
        for (u64 j = 0; j < mBits; j++)
        {
            tmp += mC[mBitMapper[j][i]] << (mBits - 1 - j);
        }

        x[i] = mLabelsRev[tmp];
    }

    //map_c_to_x(chat, xhat);
    for (u64 i = 0; i < mN; i++)
    {
        tmp = 0;
        for (u64 j = 0; j < mBits; j++)
        {
            tmp += chat[mBitMapper[j][i]] << (mBits - 1 - j);
        }

        xhat[i] = mLabelsRev[tmp];
    }

    double cw_dis_euc = 0;
    for (u64 i = 0; i < mN; i++)
    {
#ifdef ENCODE
        cw_dis_euc += (mConstellation.X()[x[i]] - mConstellation.X()[xhat[i]]) * (mConstellation.X()[x[i]] - mConstellation.X()[xhat[i]]);
#else
        cw_dis_euc += (mConstellation.X()[0] - mConstellation.X()[xhat[i]]) * (mConstellation.X()[0] - mConstellation.X()[xhat[i]]);
#endif
    }
    std::vector<u64> failed_bits_idx(cw_dis);
    j = 0;
    for (u64 i = 0; i < mLdpcCode->nc(); i++)
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
*/

} // namespace ldpc
