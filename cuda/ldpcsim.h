#pragma once

#include "ldpc/ldpc.h"

#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <chrono>
#include <vector>


#define TIME_PROF(__LOG, __EXEC, __UNIT) \
		do { \
			std::string str_unit = std::string(__UNIT);\
			float a=1;\
			if (str_unit == std::string("s")) {\
				a=1e-9;\
			} else if (str_unit == std::string("ms")) {\
				a=1e-6;\
			} else if (str_unit == std::string("us")) {\
				a=1e-3;\
			} else if (str_unit == std::string("ns")) {\
				a=1;\
			} else {\
				a=1; str_unit=std::string("ns");\
			}\
			auto start = std::chrono::high_resolution_clock::now();\
			__EXEC;\
			auto elapsed = std::chrono::high_resolution_clock::now() - start;\
			printf("[TIMEPROF]: " __LOG ": ");\
			printf("%.3f %s\n", static_cast<double>(std::chrono::duration_cast<chrono::nanoseconds>(elapsed).count())*a, str_unit.c_str());\
		} while(0);



#define MAX_FILENAME_LEN 256
#define MAX_LLR 9999.9
#define MIN_LLR -9999.9

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327
#endif

namespace ldpc
{
	class ldpc_sim
	{
	public:
		ldpc_sim(ldpc::ldpc_code* code, const char* simFileName, const char* mapFileName);
		~ldpc_sim();

		void read_bit_mapping_file(const char* filename);
		void print();
		void destroy();

		void calc_llrs(const double& y, const double& sigma2, double* llrs_out);
		double simulate_awgn(uint64_t* x, double* y, const double& sigma2);

		void log_error(bits_t* c, const uint64_t frame_num, const double snr);

		void encode() {}
		void encode_all0(uint64_t* x, bits_t* c);
		void map_c_to_x(bits_t* c, size_t* x);

		void start();

		static double randn();
	private:
		ldpc_code* mLdpcCode;
		ldpc_decoder* mLdpcDecoder;

		cstll_t* cstll;

		uint64_t n;
		uint16_t M;
		uint16_t bits;
		uint64_t max_frames;
		uint64_t min_fec;
		uint64_t bp_iter;
		double* snrs;
		uint16_t* labels;
		uint16_t* labels_rev;
		bool decoder_terminate_early;
		double SE;
		size_t num_snrs;
		char logfile[MAX_FILENAME_LEN];
		size_t** bit_mapper;
		size_t* bits_pos;
	};


	class ldpc_sim_device
	{
	public:
		ldpc_sim_device(ldpc::ldpc_code* code, const char* simFileName, const char* mapFileName);
		~ldpc_sim_device();

		void print();
		void destroy();

		void start();

		double curandn();

	private:
		thrust::device_vector<double> mSNR;
		thrust::device_vector<uint16_t> mLabels;
		thrust::device_vector<uint16_t> mLabelsRev;
	};
}
