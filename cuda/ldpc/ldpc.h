#pragma once

#include <iostream>
#include <stdint.h>

using namespace std;

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
    uint64_t n() const;
    uint64_t k() const;
    uint64_t m() const;
    uint64_t nnz() const;
    size_t *cw() const;
    size_t *vw() const;
    size_t **cn() const;
    size_t **vn() const;
    size_t *r() const;
    size_t *c() const;

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
    size_t* puncture; /* array pf punctured bit indices */
    size_t num_puncture; /* number of punctured bits */
    size_t num_puncture_sys; /* number of punctured bits in systematic part */
    size_t num_puncture_par; /* number of punctured bits in parity part */
    size_t* shorten; /* array of shortened bit indices */
    size_t num_shorten; /* number of shortened bits */
    uint64_t nct; /* number of transmitted code bits */
    uint64_t kct; /* number of transmitted information bits */
    uint64_t mct; /* number of transmitted parity check bits */
    size_t max_dc;

};
