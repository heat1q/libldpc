#pragma once

#include "../ldpc/ldpc.h"
#include "../ldpc/decoder.h"

#define MAX_FILENAME_LEN 256
#define MAX_LLR 9999.9
#define MIN_LLR -9999.9

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327
#endif

namespace pgd
{

//struct where simulation results are saved
typedef struct
{
    double* fer;
    double* ber;
    double* avg_iter;
    double* time;
    std::size_t* fec;
    std::size_t* frames;
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
    ldpc_sim(ldpc_code *pCode, const char *pSimFileName, const char *pMapFileName, std::uint16_t numThreads, std::size_t seed);
    ldpc_sim(ldpc_code *pCode, const char *pSimFileName, const char *pMapFileName, std::uint16_t numThreads, std::size_t seed, sim_results_t* results);

    void start(std::uint8_t *stopFlag);

    void simulate_awgn(double pSigma2, std::uint16_t threadid);
    void encode();
    void map();

    void print();
    void log_error(std::size_t pFrameNum, double pSNR);
    void print_file_header(const char *binaryFile, const char *codeFile, const char *simFile, const char *mapFile);

    void free_results();

    std::vector<ldpc_decoder> mLdpcDecoder;
    ldpc_code *mLdpcCode;

    const constellation &cstll() const { return mConstellation; }
    std::size_t n() const { return mN; }
    std::size_t bits() const { return mBits; }
    const vec_labels_t &labels() const { return mLabels; }
    const vec_labels_t &labels_rev() const { return mLabelsRev; }
    const vec_size_t &bits_pos() const { return mBitPos; }
    const mat_size_t &bit_mapper() const { return mBitMapper; }

private:
    constellation mConstellation;
    std::uint16_t mThreads;

    std::size_t mN;
    std::size_t mBits;
    std::size_t mMaxFrames;
    std::size_t mMinFec;
    std::size_t mBPIter;

    std::string mLogfile;
    double mSE;

    vec_double_t mSnrs;
    vec_labels_t mLabels;
    vec_labels_t mLabelsRev;
    vec_size_t mBitPos;
    mat_size_t mBitMapper;

    mat_double_t mX;
    mat_double_t mY;
    mat_bits_t mC;

    std::size_t mRNGSeed;
    std::vector<std::mt19937> mRNG;
    std::vector<std::normal_distribution<double>> mRandNormal;

    sim_results_t* mResults;
};
} // namespace pgd
