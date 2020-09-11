#pragma once

#include "channel.h"

#define MAX_FILENAME_LEN 256
#define MAX_LLR 9999.9
#define MIN_LLR -9999.9

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327
#endif

namespace ldpc
{
    enum channel_type
    {
        AWGN = 1,
        BSC,
        BEC
    };

    //struct where simulation results are saved
    typedef struct
    {
        double *fer;
        double *ber;
        double *avg_iter;
        double *time;
        u64 *fec;
        u64 *frames;
    } sim_results_t;

    class ldpc_sim
    {
    public:
        ldpc_sim() = default;
        ldpc_sim(const ldpc_code *code,
                 const param_map &decoderParams,
                 const param_map &channelParams,
                 const param_map &simulationParams);
        ldpc_sim(const ldpc_code *code,
                 const param_map &decoderParams,
                 const param_map &channelParams,
                 const param_map &simulationParams,
                 sim_results_t *results);

        void start(bool *stopFlag);

        //void log_error(u64 pFrameNum, double pSNR);
        //void print_file_header(const char *binaryFile, const char *codeFile, const char *simFile, const char *mapFile);

        friend std::ostream &operator<<(std::ostream &os, const ldpc_sim &sim);        

    private:
        const ldpc_code *mLdpcCode;
        std::vector<std::shared_ptr<ldpc_decoder>> mLdpcDecoder;
        std::vector<std::shared_ptr<channel>> mChannel;

        param_map mDecoderParams;
        param_map mChannelParams;
        param_map mSimulationParams;

        sim_results_t *mResults;
    };
} // namespace ldpc
