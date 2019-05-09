#pragma once

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

#define MAX_FILENAME_LEN 256
#define MAX_LLR 9999.9
#define MIN_LLR -9999.9

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327
#endif

namespace ldpc
{

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
	ldpc_sim(ldpc_code *pCode, ldpc_decoder *pDec, const char *pSimFileName, const char *pMapFileName);

	void start();

	double randn();
	double simulate_awgn(double pSigma2);
	void encode() {}
	void encode_all0();
	void map_c_to_x();
	void calc_llrs(double sigma2);

	void print();
	void log_error(std::size_t pFrameNum, double pSNR);
#ifdef LOG_TP
	std::size_t frame_const_time(double pSigma2, std::size_t pCount);
#endif

	ldpc_decoder *mLdpcDecoder;
	vec_size_t mX;
	vec_double_t mY;
	vec_bits_t mC;

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
} // namespace ldpc
