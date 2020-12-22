#include "ldpcsim.h"

#include <omp.h>

namespace ldpc
{
    ldpc_sim::ldpc_sim(const std::shared_ptr<ldpc_code> &code,
                       const decoder_param &decoderParams,
                       const channel_param &channelParams,
                       const simulation_param &simulationParams)
        : ldpc_sim(code, decoderParams, channelParams, simulationParams, nullptr)
    {
    }

    ldpc_sim::ldpc_sim(const std::shared_ptr<ldpc_code> &code,
                       const decoder_param &decoderParams,
                       const channel_param &channelParams,
                       const simulation_param &simulationParams,
                       sim_results_t *results)
        : mLdpcCode(code),
          mDecoderParams(decoderParams),
          mChannelParams(channelParams),
          mSimulationParams(simulationParams),
          mResults(results)
    {
        try
        {
            //results may vary with same seed, since some threads are executed more than others
            for (u32 i = 0; i < mSimulationParams.threads; ++i)
            {
                // initialize the correct channel
                if (mChannelParams.type == std::string("AWGN"))
                {
                    mChannel.push_back(
                        std::make_shared<channel_awgn>(
                            channel_awgn(
                                mLdpcCode,
                                mDecoderParams,
                                mChannelParams.seed + i,
                                1.
                            )
                        )
                    );
                }
                else if (mChannelParams.type == std::string("BSC"))
                {
                    mChannel.push_back(
                        std::make_shared<channel_bsc>(
                            channel_bsc(
                                mLdpcCode,
                                mDecoderParams,
                                mChannelParams.seed + i,
                                0.
                            )
                        )
                    );
                }
                else if (mChannelParams.type == std::string("BEC"))
                {
                    mChannel.push_back(
                        std::make_shared<channel_bec>(
                            channel_bec(
                                mLdpcCode,
                                mDecoderParams,
                                mChannelParams.seed + i,
                                0.
                            )
                        )
                    );
                }
                else
                {
                    throw std::runtime_error("No channel selected.");
                }
            }            
        }
        catch (std::exception &e)
        {
            std::cout << "Error: ldpc_sim::ldpc_sim() " << e.what() << std::endl;
            exit(EXIT_FAILURE);
        }
    }

    std::ostream &operator<<(std::ostream &os, const ldpc_sim &sim)
    {
        os << "== Decoder Parameters\n";
        os << sim.mDecoderParams << "\n";

        os << "== Channel Parameters\n";
        os << sim.mChannelParams << "\n";

        os << "== Simulation Parameters\n";
        os << sim.mSimulationParams << "\n";
        return os;
    }

    void ldpc_sim::start(bool *stopFlag)
    {
        u64 frames;
        u64 bec = 0;
        u64 fec = 0;
        u64 iters;
        u64 bec_tmp;

        ldpc::vec_double_t xVals;
        double val = mChannelParams.xRange[0];
        while (val < mChannelParams.xRange[1])
        {
            xVals.push_back(val);
            val += mChannelParams.xRange[2];
        }

        auto minFec = mSimulationParams.fec;
        auto maxFrames = mSimulationParams.maxFrames;
        std::string xValType = "SNR";

        if (mChannelParams.type == std::string("BSC") || mChannelParams.type == std::string("BEC"))
        {
            xValType = "EPS";
            // reverse the epsilon values, since we should start at the worst
            // crossover probability increase to the best
            std::reverse(xVals.begin(), xVals.end());
        }

        std::vector<std::string> printResStr(xVals.size() + 1, std::string());
        std::ofstream fp;
        char resStr[128];

        #ifndef LIB_SHARED
        #ifdef LOG_FRAME_TIME
        printResStr[0].assign("snr fer ber frames avg_iter frame_time");
        #else
        printResStr[0].assign("snr fer ber frames avg_iter");
        #endif
        #endif


        std::cout << "=============================" <<           "===========================================================" << std::endl;
        std::cout << "  FEC   |      FRAME     |   " <<  xValType << "   |    BER     |    FER     | AVGITERS  |  TIME/FRAME   \n";
        std::cout << "========+================+===" <<           "======+============+============+===========+==============" << std::endl;

        for (u64 i = 0; i < xVals.size(); ++i)
        {
            bec = 0;
            fec = 0;
            frames = 0;
            iters = 0;

            auto timeStart = std::chrono::high_resolution_clock::now();

            #pragma omp parallel default(none)                                                                        \
                num_threads(mSimulationParams.threads) private(bec_tmp)                                               \
                firstprivate(xVals, mLdpcCode, stdout)                                                                \
                shared(stopFlag, timeStart, mChannel, fec, bec, frames, printResStr, fp, resStr, i, minFec, maxFrames)\
                reduction(+: iters)
            {
                unsigned tid = omp_get_thread_num();

                // reconfigure channel to match parameter
                mChannel[tid]->set_channel_param(xVals[i]);

                do
                {
                    if (!mLdpcCode->G().empty())
                    {
                        mChannel[tid]->encode_and_map();
                    }

                    // calculate the channel transitions
                    mChannel[tid]->simulate();

                    // calculate the corresponding LLRs, depending on the channel
                    mChannel[tid]->calculate_llrs();

                    //decode
                    iters += mChannel[tid]->decode();

                    if (fec < minFec)
                    {
                        #pragma omp atomic update
                        ++frames;

                        // count the bit errors
                        bec_tmp = 0;
                        for (auto ci : mLdpcCode->bit_pos())
                        {
                            bec_tmp += (mChannel[tid]->estimate()[ci] != mChannel[tid]->codeword()[ci]);
                        }

                        if (bec_tmp > 0)
                        {
                            auto timeNow = std::chrono::high_resolution_clock::now();
                            auto timeFrame = timeNow - timeStart; //eliminate const time for printing etc
                            u64 tFrame = static_cast<u64>(std::chrono::duration_cast<std::chrono::microseconds>(timeFrame).count());
                            tFrame = tFrame / frames;
                            #pragma omp critical
                            {
                                bec += bec_tmp;
                                ++fec;

                                #ifndef LIB_SHARED
                                printf("\r %2lu/%2lu  |  %12lu  |  %.3f  |  %.2e  |  %.2e  |  %.1e  |  %.3fms",
                                       fec, minFec, frames, xVals[i],
                                       static_cast<double>(bec) / (frames * mLdpcCode->nc()), //ber
                                       static_cast<double>(fec) / frames,                     //fer
                                       static_cast<double>(iters) / frames,                   //avg iters
                                       static_cast<double>(tFrame) * 1e-3);                   //frame time tFrame
                                fflush(stdout);

                                #ifdef LOG_FRAME_TIME
                                sprintf(resStr, "%lf %.3e %.3e %lu %.3e %.6f",
                                        xVals[i], static_cast<double>(fec) / frames, static_cast<double>(bec) / (frames * mLdpcCode->nc()),
                                        frames, static_cast<double>(iters) / frames, static_cast<double>(tFrame) * 1e-6);
                                #else
                                sprintf(resStr, "%lf %.3e %.3e %lu %.3e",
                                        xVals[i], static_cast<double>(fec) / frames, static_cast<double>(bec) / (frames * mLdpcCode->nc()),
                                        frames, static_cast<double>(iters) / frames);
                                #endif
                                printResStr[i + 1].assign(resStr);

                                try
                                {
                                    fp.open(mSimulationParams.resultFile);
                                    for (const auto &x : printResStr)
                                    {
                                        fp << x << "\n";
                                    }
                                    fp.close();
                                }
                                catch (...)
                                {
                                    printf("Warning: can not open logfile for writing\n");
                                }

                                #ifdef LOG_CW
                                //log_error(frames, xVals[i]);
                                #endif
                                #endif

                                //save to result struct
                                if (mResults != nullptr)
                                {
                                    mResults->fer[i] = static_cast<double>(fec) / frames;
                                    mResults->ber[i] = static_cast<double>(bec) / (frames * mLdpcCode->nc());
                                    mResults->avg_iter[i] = static_cast<double>(iters) / frames;
                                    mResults->time[i] = static_cast<double>(tFrame) * 1e-6;
                                    mResults->fec[i] = fec;
                                    mResults->frames[i] = frames;
                                }

                                timeStart += std::chrono::high_resolution_clock::now() - timeNow; //dont measure time for printing files
                            }
                        }
                    }
                } while (fec < minFec && frames < maxFrames && !*stopFlag); //end while
            }
            #ifndef LIB_SHARED
            printf("\n");
            #endif
        } //end for

        *resStr = 0;
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
