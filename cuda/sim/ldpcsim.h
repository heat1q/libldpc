#pragma once

#include <chrono>
#include <vector>
#include <string>
#include <fstream>

#include "curand_kernel.h"

#ifndef GPU_ID
#define GPU_ID 0
#endif

#include "../ldpc/ldpc.h"


#define TIME_PROF(log, exec, unit)                                                                                                             \
	do                                                                                                                                         \
	{                                                                                                                                          \
		std::string str_unit = std::string(unit);                                                                                              \
		float a = 1;                                                                                                                           \
		if (str_unit == std::string("s"))                                                                                                      \
		{                                                                                                                                      \
			a = 1e-9;                                                                                                                          \
		}                                                                                                                                      \
		else if (str_unit == std::string("ms"))                                                                                                \
		{                                                                                                                                      \
			a = 1e-6;                                                                                                                          \
		}                                                                                                                                      \
		else if (str_unit == std::string("us"))                                                                                                \
		{                                                                                                                                      \
			a = 1e-3;                                                                                                                          \
		}                                                                                                                                      \
		else if (str_unit == std::string("ns"))                                                                                                \
		{                                                                                                                                      \
			a = 1;                                                                                                                             \
		}                                                                                                                                      \
		else                                                                                                                                   \
		{                                                                                                                                      \
			a = 1;                                                                                                                             \
			str_unit = std::string("ns");                                                                                                      \
		}                                                                                                                                      \
		auto start = std::chrono::high_resolution_clock::now();                                                                                \
		exec;                                                                                                                                  \
		auto elapsed = std::chrono::high_resolution_clock::now() - start;                                                                      \
		printf("[TIMEPROF]: " log ": ");                                                                                                       \
		printf("%.3f %s\n", static_cast<double>(std::chrono::duration_cast<std::chrono::nanoseconds>(elapsed).count()) * a, str_unit.c_str()); \
	} while (0);

#ifndef NUMK_THREADS
#define NUMK_THREADS 128
#endif

#define gpuErrchk(ans)                              \
	{                                               \
		ldpc::gpuAssert((ans), __FILE__, __LINE__); \
	}

#define MAX_FILENAME_LEN 256
#define MAX_LLR 9999.9
#define MIN_LLR -9999.9

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327
#endif

//These are used for variable sized temporary arrays inside kernel functions
#ifndef SIM_NUM_BITS
#define SIM_NUM_BITS 32
#endif
#ifndef DEC_MAX_DC
#define DEC_MAX_DC 64
#endif

namespace ldpc
{
using vec_ldpc_dec_t = cuda_vector<cuda_ptr<ldpc_decoder_device>>;
using mat_curandState_t = cuda_vector<cuda_vector<curandState_t>>;

class constellation
{
  public:
	__host__ __device__ constellation() {}
	__host__ constellation(labels_t pM);
	__host__ void mem_prefetch()
	{
		mPX.mem_prefetch();
		mX.mem_prefetch();
	}

	__host__ __device__ const vec_double_t &pX() const { return mPX; }
	__host__ __device__ const vec_double_t &X() const { return mX; }
	__host__ __device__ labels_t M() const { return mM; }
	__host__ __device__ labels_t log2M() const { return mLog2M; }

  private:
	vec_double_t mPX;
	vec_double_t mX;
	labels_t mM;
	labels_t mLog2M;
};

class ldpc_sim_device
{
  public:
	__host__ ldpc_sim_device(cuda_ptr<ldpc_code_device> &pCode, const char *pSimFileName, const char *pMapFileName, labels_t pNumThreads);

	__host__ void start_device();
	__host__ void start();

	__host__ double randn();
	__host__ __device__ double simulate_awgn(double pSigma2);
	__host__ __device__ void encode() {}
	__host__ __device__ void encode_all0();
	__host__ __device__ void map_c_to_x();
	__host__ __device__ void calc_llrs(double sigma2);

	__host__ void print();
	__host__ void log_error(std::size_t pFrameNum, double pSNR, labels_t pThreads);
#ifdef LOG_TP
	__host__ std::size_t frame_const_time(double pSigma2, std::size_t pCount);
#endif

	//changed with each frame
	vec_ldpc_dec_t mLdpcDecoderVec;
	mat_size_t mX;
	mat_double_t mY;
	mat_bits_t mC;

	mat_curandState_t mCurandState;
	mat_curandState_t mCurandStateEncoding;

	cuda_ptr<ldpc_code_device> mLdpcCode;

	__host__ __device__ const constellation &cstll() const { return mConstellation; }
	__host__ __device__ std::size_t n() const { return mN; }
	__host__ __device__ std::size_t bits() const { return mBits; }
	__host__ __device__ const vec_labels_t &labels() const { return mLabels; }
	__host__ __device__ const vec_labels_t &labels_rev() const { return mLabelsRev; }
	__host__ __device__ const vec_size_t &bits_pos() const { return mBitPos; }
	__host__ __device__ const mat_size_t &bit_mapper() const { return mBitMapper; }

  private:
	constellation mConstellation;
	labels_t mThreads;

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
};

inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort = true)
{
	if (code != cudaSuccess)
	{
		fprintf(stderr, "GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
		if (abort)
			exit(code);
	}
}
} // namespace ldpc
