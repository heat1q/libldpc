#pragma once

#include "../core/ldpc.h"
#include "../decoding/decoder.h"

#define MAX_FILENAME_LEN 256
#define MAX_LLR 9999.9
#define MIN_LLR -9999.9

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327
#endif

namespace ldpc
{

enum channel_type {
    AWGN = 1,
    BSC,
    BEC
};

//struct where simulation results are saved
typedef struct
{
    double* fer;
    double* ber;
    double* avg_iter;
    double* time;
    u64* fec;
    u64* frames;
} sim_results_t;


class ldpc_sim
{
public:
    ldpc_sim() = default;
    ldpc_sim(const ldpc_code *code,
             const std::string &outFile,
             const vec_double_t &snrRange,
             const unsigned numThreads,
             const u64 seed,
             const enum channel_type channel,
             const unsigned iters,
             const u64 maxFrames,
             const u64 fec,
             const bool earlyTerm);
    ldpc_sim(const ldpc_code *code,
             const std::string &outFile,
             const vec_double_t &snrRange,
             const unsigned numThreads,
             const u64 seed,
             const enum channel_type channel,
             const unsigned iters,
             const u64 maxFrames,
             const u64 fec,
             const bool earlyTerm,
             sim_results_t *results);

    void start(bool *stopFlag);

    void simulate_awgn(double pSigma2, unsigned threadid);

    //void log_error(u64 pFrameNum, double pSNR);
    //void print_file_header(const char *binaryFile, const char *codeFile, const char *simFile, const char *mapFile);

    friend std::ostream &operator<<(std::ostream &os, const ldpc_sim &sim);

    const vec_u64& bit_pos() const { return mBitPos; }
    
private:
    const ldpc_code *mLdpcCode;
    const std::string mLogfile;
    std::vector<ldpc_decoder> mLdpcDecoder;

    const unsigned mThreads;

    const unsigned mBPIter;
    const u64 mMaxFrames;
    const u64 mMinFec;

    vec_u64 mBitPos; // position of transmitted bits, i.e. puncture/shorten exluded
    vec_double_t mSnrs;
    const enum channel_type mChannel;

    //channel i/o
    mat_double_t mX;
    mat_double_t mY;
    mat_bits_t mC;

    const u64 mRNGSeed;
    std::vector<std::mt19937> mRNG;
    std::vector<std::normal_distribution<double>> mRandNormal;

    sim_results_t* mResults;
};
} // namespace ldpc
