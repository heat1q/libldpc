#pragma once

#define QC_LYR_DEC

#include <iostream>
#include <stdint.h>

typedef uint8_t bits_t;
typedef uint16_t labels_t;
typedef uint32_t symbols_t;

namespace ldpc {

class Ldpc_Code_cl
{
public:
    Ldpc_Code_cl(const char* filename);
    virtual ~Ldpc_Code_cl();

    #ifdef QC_LYR_DEC
	Ldpc_Code_cl();
    Ldpc_Code_cl(const char *filename, const char* clfile);
    void setup_layers(const char* clfile);
	void setup_layers_managed(const char* clfile);
	void setup_code_managed(const char *filename, const char* clfile);
	void destroy_ldpc_code_managed();
    #endif

    void print_ldpc_code();

    void destroy_ldpc_code();


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

    #ifdef QC_LYR_DEC
    __host__ __device__ uint64_t nl() const;
    __host__ __device__ uint64_t* lw() const;
    __host__ __device__ uint64_t** layers() const;
    #endif

private:
	uint8_t level = 0;
    uint64_t n_c;
    uint64_t k_c;
    uint64_t m_c;
    uint64_t nnz_c;
    size_t* cw_c; /* denotes the check weight of each check node, i.e., # of connected VN; dimensions cw[mc] */
    size_t* vw_c; /* denotes the variable weight, i.e., # of connected CN; dimensions vw[nc] */
    size_t** cn_c; /* denotes the check neighbors, i.e. connected VN, for each check node as index in c/r; dimensions cn[mc][cw[i]] */
    size_t** vn_c; /* denotes the var neighbors, i.e., connected CN, for each variable node as index in c/r; dimensions vn[nc][vw[i]] */
    size_t* r_c; /* non zero row indices; length nnz */
    size_t* c_c; /* non zero check indices; length nnz */
    size_t* puncture_c; /* array pf punctured bit indices */
    size_t num_puncture_c; /* number of punctured bits */
    size_t num_puncture_sys_c; /* number of punctured bits in systematic part */
    size_t num_puncture_par_c; /* number of punctured bits in parity part */
    size_t* shorten_c; /* array of shortened bit indices */
    size_t num_shorten_c; /* number of shortened bits */
    uint64_t nct_c; /* number of transmitted code bits */
    uint64_t kct_c; /* number of transmitted information bits */
    uint64_t mct_c; /* number of transmitted parity check bits */
    size_t max_dc_c;
    #ifdef QC_LYR_DEC
    uint64_t nl_c; //number of layers
    uint64_t* lw_c; //layer weight
    uint64_t** layers_c;
    #endif
};


class Ldpc_Decoder_cl
{
public:
	__host__ __device__ Ldpc_Decoder_cl();
	Ldpc_Decoder_cl(Ldpc_Code_cl* code);
	__host__ __device__ ~Ldpc_Decoder_cl();

	void setup_decoder(Ldpc_Code_cl* code);
	__device__ void setup_decoder_device(Ldpc_Code_cl* code);

	__host__ __device__ void destroy_dec();
    __host__ __device__ bool is_codeword();

    uint64_t decode_legacy(double* llr_in, double* llr_out, const uint64_t& max_iter, const bool& early_termination);
    uint64_t decode_layered_legacy(double* llr_in, double* llr_out, const uint64_t& max_iter, const bool& early_termination);

    Ldpc_Code_cl* ldpc_code;
    double* l_v2c;
    double* l_c2v;
    double* f;
    double* b;
    double* lsum;
    double* l_c2v_pre;

    bits_t* c_out;

private:
    	bool init = false;
};


template<typename T>
extern void printVector(T* x, const size_t& l);

void dec2bin(uint64_t val, uint8_t m);
__host__ __device__ double jacobian(const double& L1, const double& L2);
__host__ __device__ double jacobian_lin_approx(const double& L);
__host__ __device__ int8_t sign(const double& a);


namespace cudakernel {
	__global__ void setup_decoder(Ldpc_Code_cl* code_managed, Ldpc_Decoder_cl** dec_ptr);
	__global__ void destroy_decoder(Ldpc_Decoder_cl** dec_ptr);
    __global__ void clean_decoder(Ldpc_Decoder_cl** dec_ptr);
    __global__ void decode_lyr_vnupdate(Ldpc_Decoder_cl** dec_ptr, double* llr_in, size_t i_nnz);
    __global__ void decode_lyr_cnupdate(Ldpc_Decoder_cl** dec_ptr, size_t i_nnz, uint64_t L);
    __global__ void decode_lyr_sumllr(Ldpc_Decoder_cl** dec_ptr, size_t i_nnz);
    __global__ void decode_lyr_appcalc(Ldpc_Decoder_cl** dec_ptr, double* llr_in, double* llr_out);
}

}
