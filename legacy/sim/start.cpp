#include "ldpcsim.h"

namespace ldpc
{
//start simulation on cpu
void ldpc_sim::start()
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
#ifdef LOG_TP
        std::size_t tconst = frame_const_time(sigma2, mMinFec);
#endif
        auto timeStart = std::chrono::high_resolution_clock::now();
        do
        {
            encode_all0();
            simulate_awgn(sigma2);

            //puncturing & shortening
            if (mLdpcCode->num_puncture() != 0)
            {
                for (std::size_t j = 0; j < mLdpcCode->num_puncture(); j++)
                {
                    mLdpcDecoder->mLLRIn[mLdpcCode->puncture()[j]] = 0;
                }
            }
            if (mLdpcCode->num_shorten() != 0)
            {
                for (std::size_t j = 0; j < mLdpcCode->num_shorten(); j++)
                {
                    mLdpcDecoder->mLLRIn[mLdpcCode->shorten()[j]] = 99999.9;
                }
            }

            calc_llrs(sigma2);
            //decode
            iters += mLdpcDecoder->decode_layered();

            ++frames;

            bec_tmp = 0;
            for (std::size_t j = 0; j < mLdpcCode->nc(); ++j)
            {
                bec_tmp += (mLdpcDecoder->mLLROut[j] <= 0);
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
                std::size_t tDec = std::max(1, static_cast<int>(tFrame) - static_cast<int>(tconst));
                printf("\r %2lu/%2lu  |  %12lu  |  %.3f  |  %.2e  |  %.2e  |  %.1e  |  %7.3fms  |  %6luus  |  %.2fMbits/s",
                       fec, mMinFec, frames, mSnrs[i],
                       static_cast<double>(bec) / (frames * mLdpcCode->nc()), //ber
                       static_cast<double>(fec) / frames,                     //fer
                       static_cast<double>(iters) / frames,                   //avg iters
                       static_cast<double>(tFrame) * 1e-3,                    //frame time
                       tDec,                                                  //decoding time
                       static_cast<double>(mLdpcCode->nc()) / (tDec));        //decoding throughput
#else
                printf("\r %2lu/%2lu  |  %12lu  |  %.3f  |  %.2e  |  %.2e  |  %.1e  |  %.3fms",
                       fec, mMinFec, frames, mSnrs[i],
                       static_cast<double>(bec) / (frames * mLdpcCode->nc()), //ber
                       static_cast<double>(fec) / frames,                     //fer
                       static_cast<double>(iters) / frames,                   //avg iters
                       static_cast<double>(tFrame) * 1e-3);                   //frame time
#endif
                fflush(stdout);

#ifdef LOG_FRAME_TIME
                sprintf(resStr, "%lf %.3e %.3e %lu %.3e %.6f",
                        mSnrs[i], static_cast<double>(fec) / frames, static_cast<double>(bec) / (frames * mLdpcCode->nc()),
                        frames, static_cast<double>(iters) / frames, static_cast<double>(tFrame) * 1e-6);
#elif defined LOG_TP
                sprintf(resStr, "%lf %.3e %.3e %lu %.3e %.6f %.6f %lu",
                        mSnrs[i], static_cast<double>(fec) / frames, static_cast<double>(bec) / (frames * mLdpcCode->nc()),
                        frames, static_cast<double>(iters) / frames, static_cast<double>(tFrame) * 1e-6,
                        static_cast<double>(tDec) * 1e-6, static_cast<std::size_t>(mLdpcCode->nc() / ((tFrame - tconst) * 1e-6)));
#else
                sprintf(resStr, "%lf %.3e %.3e %lu %.3e",
                        mSnrs[i], static_cast<double>(fec) / frames, static_cast<double>(bec) / (frames * mLdpcCode->nc()),
                        frames, static_cast<double>(iters) / frames);
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

                log_error(frames, mSnrs[i]);

                timeStart += std::chrono::high_resolution_clock::now() - timeNow; //dont measure time for printing files
            }
        } while (fec < mMinFec && frames < mMaxFrames); //end while
        printf("\n");
    } //end for
}
} // namespace ldpc