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
    Ldpc_Code_cl(const char *filename, const char* clfile);
    void setup_layers(const char* clfile);
    #endif

    void print_ldpc_code();

    void destroy_ldpc_code();


    //getter functions
    uint64_t nc() const;
    uint64_t kc() const;
    uint64_t mc() const;
    uint64_t nnz() const;
    size_t *cw() const;
    size_t *vw() const;
    size_t **cn() const;
    size_t **vn() const;
    size_t *r() const;
    size_t *c() const;

    uint64_t nct() const;
    uint64_t kct() const;
    uint64_t mct() const;

    size_t *puncture() const;
    size_t num_puncture() const;
    size_t *shorten() const;
    size_t num_shorten() const;
    size_t max_dc() const;

    #ifdef QC_LYR_DEC
    uint64_t nl() const;
    uint64_t* lw() const;
    uint64_t** layers() const;
    #endif

private:
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
	Ldpc_Decoder_cl();
	Ldpc_Decoder_cl(Ldpc_Code_cl* code);
	virtual ~Ldpc_Decoder_cl();

	void setup_decoder(Ldpc_Code_cl* code);

	void destroy_dec();
    bool is_codeword_global(bits_t *c);

    uint64_t decode(double* llr_in, double* llr_out, const uint64_t& max_iter, const bool& early_termination);
    #ifdef QC_LYR_DEC
    uint64_t decode_layered(double* llr_in, double* llr_out, const uint64_t& MaxIter, const bool& early_termination);
    void decode_lyr_nodeupdate_global(double* llr_in);
    void decode_lyr_sumllr_global();
    void decode_lyr_appcalc_global(double* llr_in, double* llr_out);
    #endif
private:
    Ldpc_Code_cl* ldpc_code;
    double* l_v2c;
    double* l_c2v;
    double* f;
    double* b;
    double* lsum;

    bits_t* c_out;
    bits_t* synd;
};


template<typename T>
extern void printVector(T* x, const size_t& l);

void dec2bin(uint64_t val, uint8_t m);
double jacobian(const double& L1, const double& L2);
double jacobian_lin_approx(const double& L);
int8_t sign(const double& a);
}
