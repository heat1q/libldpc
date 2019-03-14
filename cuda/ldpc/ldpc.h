#pragma once

#include <iostream>
#include <stdint.h>

using namespace std;

typedef uint8_t bits_t;
typedef uint16_t labels_t;
typedef uint32_t symbols_t;

class LDPC
{
public:
    LDPC(const char *filename);
    ~LDPC();

    bool read_file(const char *filename);
    void print_ldpc_code();

    void encode_all0();
    void decode_layered();





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

};
