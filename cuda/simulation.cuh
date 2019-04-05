#pragma once

#include <chrono>
#include "ldpc/ldpc.cuh"


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

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
	if (code != cudaSuccess)
	{
		fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
		if (abort) exit(code);
	}
}

#define MAX_FILENAME_LEN 256
#define MAX_LLR 9999.9
#define MIN_LLR -9999.9

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327
#endif

struct {
    double *pX;
    double *X;
    double *A;
    uint16_t M;
    uint16_t log2M;
} typedef cstll_t;


class Sim_AWGN_cl
{
public:
    Sim_AWGN_cl(ldpc::Ldpc_Code_cl* code, const char* simFileName, const char* mapFileName);
    ~Sim_AWGN_cl();

    void read_bit_mapping_file(const char* filename);
    void print();
    void destroy();

    void calc_llrs(const double& y, const double& sigma2, double* llrs_out);
    double simulate_awgn(uint64_t* x, double* y, const double& sigma2);
    __host__ __device__ static double randn();

    void encode() {}
    void encode_all0(uint64_t* x, bits_t* c);
    void map_c_to_x(bits_t* c, size_t* x);

    void start();

private:
    ldpc::Ldpc_Code_cl* ldpc_code;
	ldpc::Ldpc_Decoder_cl* ldpc_dec;

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
