#pragma once

#include "../device/cudamgd.h"

namespace ldpc
{
	typedef uint8_t bits_t;
	typedef uint16_t labels_t;
	typedef uint32_t symbols_t;

	struct {
	    double* pX;
	    double* X;
	    double* A;
	    uint16_t M;
	    uint16_t log2M;
	} typedef cstll_t;

	class ldpc_code : public cuda_mgd
	{
	public:
		ldpc_code(const char* pFileName, const char* pClFile, const bool pMgd);
		virtual ~ldpc_code();

		void setup(const char* pFileName);
		void setup_layers(const char* pClFile);

		void setup_mgd(const char* pFileName);
		void setup_layers_mgd(const char* pClFile);
		void prefetch();

		void destroy();
		void destroy_mgd();

		void print();

		//getter functions
		__host__ __device__ uint64_t nc() const;
		__host__ __device__ uint64_t kc() const;
		__host__ __device__ uint64_t mc() const;
		__host__ __device__ uint64_t nnz() const;
		__host__ __device__ size_t *cw() const;
		__host__ __device__ size_t *vw() const;
		__host__ __device__ size_t **cn() const;
		__host__ __device__ size_t **vn() const;
		__host__ __device__ size_t *r() const;
		__host__ __device__ size_t *c() const;
		__host__ __device__ uint64_t nct() const;
		__host__ __device__ uint64_t kct() const;
		__host__ __device__ uint64_t mct() const;
		__host__ __device__ size_t *puncture() const;
		__host__ __device__ size_t num_puncture() const;
		__host__ __device__ size_t *shorten() const;
		__host__ __device__ size_t num_shorten() const;
		__host__ __device__ size_t max_dc() const;
		__host__ __device__ uint64_t nl() const;
		__host__ __device__ uint64_t* lw() const;
		__host__ __device__ uint64_t** layers() const;

	private:
		uint64_t mN;
		uint64_t mK;
		uint64_t mM;
		uint64_t mNNZ;
		size_t* mCW; /* denotes the check weight of each check node, i.e., # of connected VN; dimensions cw[mc] */
		size_t* mVW; /* denotes the variable weight, i.e., # of connected CN; dimensions vw[nc] */
		size_t** mCN; /* denotes the check neighbors, i.e. connected VN, for each check node as index in c/r; dimensions cn[mc][cw[i]] */
		size_t** mVN; /* denotes the var neighbors, i.e., connected CN, for each variable node as index in c/r; dimensions vn[nc][vw[i]] */
		size_t* mR; /* non zero row indices; length nnz */
		size_t* mC; /* non zero check indices; length nnz */
		size_t* mPuncture; /* array pf punctured bit indices */
		size_t mNumPuncture; /* number of punctured bits */
		size_t mNumPunctureSys; /* number of punctured bits in systematic part */
		size_t mNumPuncturePar; /* number of punctured bits in parity part */
		size_t* mShorten; /* array of shortened bit indices */
		size_t mNumShorten; /* number of shortened bits */
		uint64_t mNCT; /* number of transmitted code bits */
		uint64_t mKCT; /* number of transmitted information bits */
		uint64_t mMCT; /* number of transmitted parity check bits */
		size_t mMaxDC;
		uint64_t mNL; //number of layers
		uint64_t* mLW; //layer weight
		uint64_t** mLayers;
	};


	class ldpc_decoder : public cuda_mgd
	{
	public:
		ldpc_decoder(ldpc_code* pCode, const uint16_t pI, const bool pEarlyTerm);
		~ldpc_decoder();

		void setup();
		void prefetch();
		void destroy();

		__host__ __device__ bool is_codeword();
		__host__ __device__ bool is_codeword_legacy();

		uint16_t decode_legacy();
		uint16_t decode_layered_legacy();
		uint16_t decode_layered();

		__host__ __device__ uint16_t max_iter() const;
		__host__ __device__ bool early_termination() const;

		ldpc_code* mLdpcCode;

		double* mLv2c;
		double* mLc2v;
		double* mF;
		double* mB;
		double* mLSum;
		double* mLc2vPre;

		double* mLLRIn;
		double* mLLROut;

		bits_t* mSynd;
		bits_t* mCO;

		uint16_t mIter;
		void* FBREF;

		bool mIsCW;
	private:
		uint16_t mMaxIter;
		bool mEarlyTerm;
	};


	//template<typename T>
	//extern void printVector(T* x, const size_t& l);

	void dec2bin(uint64_t val, uint8_t m);
	__host__ __device__ double jacobian(const double& L1, const double& L2);
	__host__ __device__ double jacobian_lin_approx(const double& L);
	__host__ __device__ int8_t sign(const double& a);


	__host__ __device__ inline const uint64_t get_num_size(const uint64_t length, const uint16_t blockSize) { return ceil((length+blockSize-1)/blockSize); }
	__host__ __device__ inline const dim3& get_gridsize_2d(const dim3& length, const dim3& blockSize) { return dim3(ceil((length.x+blockSize.x-1)/blockSize.x), ceil((length.y+blockSize.y-1)/blockSize.y)); }


	/*
	* Kernels
	*/

	namespace cudakernel
	{
		namespace decoder
		{
			__global__ void clean_decoder(ldpc_decoder* pDecMgd);
			__global__ void decode_layered(ldpc_decoder* pDecMgd);
			__global__ void decode_lyr_vnupdate(ldpc_decoder* pDecMgd, size_t pI);
			__global__ void decode_lyr_cnupdate(ldpc_decoder* pDecMgd, size_t pI, uint64_t pL);
			__global__ void decode_lyr_sumllr(ldpc_decoder* pDecMgd, size_t pI);
			__global__ void decode_lyr_appcalc(ldpc_decoder* pDecMgd);
			__global__ void calc_synd(ldpc_decoder* pDecMgd);
		}
	}

}
/*
template<typename T> void ldpc::printVector(T *x, const size_t &l)
{
    std::cout << "[";
    for (size_t i = 0; i < l-1; ++i)
        std::cout << x[i] << " ";
    std::cout << x[l-1] << "]";
}
*/
