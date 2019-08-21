#pragma once

#include "../device/cuda_container.h"

namespace ldpc
{
using bits_t = unsigned char;
using labels_t = unsigned short;
using symbols_t = unsigned;

using vec_bits_t = cuda_vector<bits_t>;
using vec_labels_t = cuda_vector<unsigned short>;
using vec_symbols_t = cuda_vector<unsigned>;
using vec_size_t = cuda_vector<std::size_t>;
using vec_double_t = cuda_vector<double>;

using mat_bits_t = cuda_vector<cuda_vector<bits_t>>;
using mat_size_t = cuda_vector<cuda_vector<std::size_t>>;
using mat_double_t = cuda_vector<cuda_vector<double>>;

class ldpc_code_device
{
public:
    __host__ ldpc_code_device(const char *pFileName, const char *pClFile);
    __host__ void mem_prefetch();
    __host__ void print();

    //getter functions
    __host__ __device__ std::size_t nc() const { return mN; };
    __host__ __device__ std::size_t kc() const { return mK; };
    __host__ __device__ std::size_t mc() const { return mM; };
    __host__ __device__ std::size_t nnz() const { return mNNZ; };
    __host__ __device__ const vec_size_t &cw() const { return mCW; };
    __host__ __device__ const vec_size_t &vw() const { return mVW; };
    __host__ __device__ const mat_size_t &cn() const { return mCN; };
    __host__ __device__ const mat_size_t &vn() const { return mVN; };
    __host__ __device__ const vec_size_t &r() const { return mR; };
    __host__ __device__ const vec_size_t &c() const { return mC; };
    __host__ __device__ std::size_t nct() const { return mNCT; };
    __host__ __device__ std::size_t kct() const { return mKCT; };
    __host__ __device__ std::size_t mct() const { return mMCT; };
    __host__ __device__ const vec_size_t &puncture() const { return mPuncture; };
    __host__ __device__ std::size_t num_puncture() const { return mNumPuncture; };
    __host__ __device__ const vec_size_t &shorten() const { return mShorten; };
    __host__ __device__ std::size_t num_shorten() const { return mNumShorten; };
    __host__ __device__ std::size_t max_dc() const { return mMaxDC; };
    __host__ __device__ std::size_t nl() const { return mNL; };
    __host__ __device__ const vec_size_t &lw() const { return mLW; };
    __host__ __device__ const mat_size_t &layers() const { return mLayers; };

private:
    std::size_t mN;
    std::size_t mK;
    std::size_t mM;
    std::size_t mNNZ;
    vec_size_t mCW;              /* denotes the check weight of each check node, i.e., # of connected VN; dimensions cw[mc] */
    vec_size_t mVW;              /* denotes the variable weight, i.e., # of connected CN; dimensions vw[nc] */
    mat_size_t mCN;              /* denotes the check neighbors, i.e. connected VN, for each check node as index in c/r; dimensions cn[mc][cw[i]] */
    mat_size_t mVN;              /* denotes the var neighbors, i.e., connected CN, for each variable node as index in c/r; dimensions vn[nc][vw[i]] */
    vec_size_t mR;               /* non zero row indices; length nnz */
    vec_size_t mC;               /* non zero check indices; length nnz */
    vec_size_t mPuncture;        /* array pf punctured bit indices */
    std::size_t mNumPuncture;    /* number of punctured bits */
    std::size_t mNumPunctureSys; /* number of punctured bits in systematic part */
    std::size_t mNumPuncturePar; /* number of punctured bits in parity part */
    vec_size_t mShorten;         /* array of shortened bit indices */
    std::size_t mNumShorten;     /* number of shortened bits */
    std::size_t mNCT;            /* number of transmitted code bits */
    std::size_t mKCT;            /* number of transmitted information bits */
    std::size_t mMCT;            /* number of transmitted parity check bits */
    std::size_t mMaxDC;
    std::size_t mNL; //number of layers
    vec_size_t mLW;  //layer weight
    mat_size_t mLayers;
};

class ldpc_decoder
{
public:
    __host__ ldpc_decoder(cuda_ptr<ldpc_code_device> &pCode, std::size_t pI, bool pEarlyTerm);
    __host__ void mem_prefetch();

    __host__ __device__ bool is_codeword();

    __host__ __device__ std::size_t max_iter() const { return mMaxIter; }
    __host__ __device__ bool early_termination() const { return mEarlyTerm; }

    cuda_ptr<ldpc_code_device> mLdpcCode;

    vec_double_t mLv2c;
    vec_double_t mLc2v;
    vec_double_t mExMsgCN;

    vec_double_t mLLRIn;
    vec_double_t mLLROut;

    vec_bits_t mSynd;
    vec_bits_t mCO;

    std::size_t mIter;

    bool mIsCW;

private:
    std::size_t mMaxIter;
    bool mEarlyTerm;
};

__host__ __device__ void dec2bin(std::size_t val, uint8_t m);
__host__ __device__ double jacobian(double L1, double L2);
__host__ __device__ double jacobian_lin_approx(double L);
__host__ __device__ int sign(double a);

__host__ __device__ inline const std::size_t get_num_size(std::size_t length, labels_t blockSize) { return ceil((length + blockSize - 1) / blockSize); }
__host__ __device__ inline const dim3 &get_gridsize_2d(const dim3 &length, const dim3 &blockSize) { return dim3(ceil((length.x + blockSize.x - 1) / blockSize.x), ceil((length.y + blockSize.y - 1) / blockSize.y)); }

/*
	* Kernels
	*/
class ldpc_sim_device;

namespace cudakernel
{
namespace decoder
{
__global__ void init_decoder(ldpc_decoder *pDecMgd);
__global__ void decode_lyr_cnupdate(ldpc_decoder *pDecMgd, std::size_t pL);
__global__ void decode_lyr_appcalc(ldpc_decoder *pDecMgd);
__global__ void calc_synd(ldpc_decoder *pDecMgd);
} // namespace decoder

namespace sim
{
__global__ void setup_rng(ldpc_sim_device *pSim);
__global__ void frame_proc(ldpc_sim_device *pSim, double pSigma2);
#ifdef LOG_TP
__global__ void frame_time(ldpc_sim_device *pSim, double pSigma2);
#endif
__global__ void encode_all0(ldpc_sim_device *pSim, labels_t pBlockID);
__global__ void awgn(ldpc_sim_device *pSim, double pSigma2, labels_t pBlockID);
__global__ void calc_llrs(ldpc_sim_device *pSim, double pSigma2, labels_t pBlockID);
__global__ void calc_llrin(ldpc_sim_device *pSim, labels_t pBlockID);
__global__ void map_c_to_x(ldpc_sim_device *pSim, labels_t pBlockID);
} // namespace sim
} // namespace cudakernel

} // namespace ldpc