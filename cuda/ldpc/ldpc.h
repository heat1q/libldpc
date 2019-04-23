#pragma once

#include "../device/cudamgd.h"
#include "../device/vectormgd.h"

namespace ldpc
{
	using bits_t = unsigned char;
	using labels_t = unsigned short;
	using symbols_t = unsigned;

	using vec_bits_t = vector_mgd<bits_t>;
	using vec_labels_t = vector_mgd<unsigned short>;
	using vec_symbols_t = vector_mgd<unsigned>;
	using vec_size_t = vector_mgd<size_t>;
	using vec_double_t = vector_mgd<double>;

	using mat_size_t = vector_mgd< vector_mgd<size_t> >;

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


	class ldpc_code_device
	{
	public:
		__host__ ldpc_code_device(const char* pFileName, const char* pClFile, bool pUseLayer);
		__host__ void mem_prefetch();
		__host__ void print();

		//getter functions
		__host__ __device__ size_t nc() const { return mN; };
		__host__ __device__ size_t kc() const { return mK; };
		__host__ __device__ size_t mc() const { return mM; };
		__host__ __device__ size_t nnz() const { return mNNZ; };
		__host__ __device__ const vec_size_t& cw() const { return mCW; };
		__host__ __device__ const vec_size_t& vw() const { return mVW; };
		__host__ __device__ const mat_size_t& cn() const { return mCN; };
		__host__ __device__ const mat_size_t& vn() const { return mVN; };
		__host__ __device__ const vec_size_t& r() const { return mR; };
		__host__ __device__ const vec_size_t& c() const { return mC; };
		__host__ __device__ size_t nct() const { return mNCT; };
		__host__ __device__ size_t kct() const { return mKCT; };
		__host__ __device__ size_t mct() const { return mMCT; };
		__host__ __device__ const vec_size_t& puncture() const { return mPuncture; };
		__host__ __device__ size_t num_puncture() const { return mNumPuncture; };
		__host__ __device__ const vec_size_t& shorten() const { return mShorten; };
		__host__ __device__ size_t num_shorten() const { return mNumShorten; };
		__host__ __device__ size_t max_dc() const { return mMaxDC; };
		__host__ __device__ size_t nl() const { return mNL; };
		__host__ __device__ const vec_size_t& lw() const { return mLW; };
		__host__ __device__ const mat_size_t& layers() const { return mLayers; };

	private:
		size_t mN;
		size_t mK;
		size_t mM;
		size_t mNNZ;
		vec_size_t mCW; /* denotes the check weight of each check node, i.e., # of connected VN; dimensions cw[mc] */
		vec_size_t mVW; /* denotes the variable weight, i.e., # of connected CN; dimensions vw[nc] */
		mat_size_t mCN; /* denotes the check neighbors, i.e. connected VN, for each check node as index in c/r; dimensions cn[mc][cw[i]] */
		mat_size_t mVN; /* denotes the var neighbors, i.e., connected CN, for each variable node as index in c/r; dimensions vn[nc][vw[i]] */
		vec_size_t mR; /* non zero row indices; length nnz */
		vec_size_t mC; /* non zero check indices; length nnz */
		vec_size_t mPuncture; /* array pf punctured bit indices */
		size_t mNumPuncture; /* number of punctured bits */
		size_t mNumPunctureSys; /* number of punctured bits in systematic part */
		size_t mNumPuncturePar; /* number of punctured bits in parity part */
		vec_size_t mShorten; /* array of shortened bit indices */
		size_t mNumShorten; /* number of shortened bits */
		size_t mNCT; /* number of transmitted code bits */
		size_t mKCT; /* number of transmitted information bits */
		size_t mMCT; /* number of transmitted parity check bits */
		size_t mMaxDC;
		size_t mNL; //number of layers
		vec_size_t mLW; //layer weight
		mat_size_t mLayers;
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

		__host__ __device__ uint16_t max_iter() const { return mMaxIter; }
		__host__ __device__ bool early_termination() const { return mEarlyTerm; }

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

	class ldpc_decoder_device
	{
	public:
		__host__ ldpc_decoder_device(cudamgd_ptr<ldpc_code_device>& pCode, size_t pI, bool pEarlyTerm);
		__host__ ldpc_decoder_device(const ldpc_decoder_device& pCopy);
		//__host__ ldpc_decoder_device(ldpc_decoder_device&& pMove) noexcept;
		__host__ ~ldpc_decoder_device();
		__host__ ldpc_decoder_device& operator=(ldpc_decoder_device pCopy) noexcept;
		__host__ void mem_prefetch();

		__host__ __device__ bool is_codeword();
		__host__ __device__ size_t decode_layered();

		__host__ __device__ size_t decode_legacy();
		__host__ __device__ bool is_codeword_legacy();

		__host__ __device__ size_t max_iter() const { return mMaxIter; }
		__host__ __device__ bool early_termination() const { return mEarlyTerm; }

		cudamgd_ptr<ldpc_code_device> mLdpcCode;

		vec_double_t mLv2c;
		vec_double_t mLc2v;
		vec_double_t mLSum;
		vec_double_t mLc2vPre;
		vec_double_t mF;
		vec_double_t mB;

		vec_double_t mLLRIn;
		vec_double_t mLLROut;

		vec_bits_t mSynd;
		vec_bits_t mCO;

		size_t mIter;
		void* FBREF;

		bool mIsCW;
	private:
		size_t mMaxIter;
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
	class ldpc_sim_device;

	namespace cudakernel
	{
		namespace decoder
		{
			__global__ void clean_decoder(ldpc_decoder_device* pDecMgd);
			__global__ void decode_layered(ldpc_decoder_device* pDecMgd);
			__global__ void decode_lyr_vnupdate(ldpc_decoder_device* pDecMgd, size_t pI);
			__global__ void decode_lyr_cnupdate(ldpc_decoder_device* pDecMgd, size_t pI, uint64_t pL);
			__global__ void decode_lyr_sumllr(ldpc_decoder_device* pDecMgd, size_t pI);
			__global__ void decode_lyr_appcalc(ldpc_decoder_device* pDecMgd);
			__global__ void calc_synd(ldpc_decoder_device* pDecMgd);
			/*
			__global__ void clean_decoder(ldpc_decoder* pDecMgd);
			__global__ void decode_layered(ldpc_decoder* pDecMgd);
			__global__ void decode_lyr_vnupdate(ldpc_decoder* pDecMgd, size_t pI);
			__global__ void decode_lyr_cnupdate(ldpc_decoder* pDecMgd, size_t pI, uint64_t pL);
			__global__ void decode_lyr_sumllr(ldpc_decoder* pDecMgd, size_t pI);
			__global__ void decode_lyr_appcalc(ldpc_decoder* pDecMgd);
			__global__ void calc_synd(ldpc_decoder* pDecMgd);
			*/
		}

		namespace sim
		{
			__global__ void setup_randn(cudamgd_ptr<ldpc_sim_device> pSim);
			__global__ void encode_all0()
			__global__ void awgn(cudamgd_ptr<ldpc_sim_device> pSim, double sigma2);
			__global__ void calc_llrs(ldpc_sim_device* pSim, double pSigma2);
			__global__ void frame_proc(cudamgd_ptr<ldpc_sim_device> pSim);
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
