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

class constellation
{
public:
    constellation() {}
    constellation(labels_t pM);

    const vec_double_t &pX() const { return mPX; }
    const vec_double_t &X() const { return mX; }
    labels_t M() const { return mM; }
    labels_t log2M() const { return mLog2M; }

private:
    vec_double_t mPX;
    vec_double_t mX;
    labels_t mM;
    labels_t mLog2M;
};

class ldpc_sim
{
public:
    ldpc_sim(ldpc_code *pCode, const char *pSimFileName, const char *pMapFileName, unsigned numThreads, u64 seed);
    ldpc_sim(ldpc_code *pCode, const char *pSimFileName, const char *pMapFileName, unsigned numThreads, u64 seed, sim_results_t* results);

    void start(bool *stopFlag);

    void simulate_awgn(double pSigma2, unsigned threadid);
    void encode();
    void map();

    void print();
    void log_error(u64 pFrameNum, double pSNR);
    void print_file_header(const char *binaryFile, const char *codeFile, const char *simFile, const char *mapFile);

    void free_results();

    std::vector<ldpc_decoder> mLdpcDecoder;
    ldpc_code *mLdpcCode;

    const constellation &cstll() const { return mConstellation; }
    u64 n() const { return mN; }
    u64 bits() const { return mBits; }
    const vec_labels_t &labels() const { return mLabels; }
    const vec_labels_t &labels_rev() const { return mLabelsRev; }
    const vec_u64 &bits_pos() const { return mBitPos; }
    const mat_u64 &bit_mapper() const { return mBitMapper; }

private:
    constellation mConstellation;
    unsigned mThreads;

    u64 mN;
    u64 mBits;
    u64 mMaxFrames;
    u64 mMinFec;
    u64 mBPIter;

    std::string mLogfile;
    double mSE;

    vec_double_t mSnrs;
    vec_labels_t mLabels;
    vec_labels_t mLabelsRev;
    vec_u64 mBitPos;
    mat_u64 mBitMapper;

    mat_double_t mX;
    mat_double_t mY;
    mat_bits_t mC;

    u64 mRNGSeed;
    std::vector<std::mt19937> mRNG;
    std::vector<std::normal_distribution<double>> mRandNormal;

    sim_results_t* mResults;
};
} // namespace ldpc
