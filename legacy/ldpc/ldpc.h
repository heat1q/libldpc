#pragma once

#include "../core/functions.h"


namespace pgd
{


/**
 * @brief LDPC code class
 * 
 */
class ldpc_code
{
public:
    ldpc_code(const char *pFileName);
    void print();

    //getter functions
    std::size_t nc() const { return mN; };
    std::size_t kc() const { return mK; };
    std::size_t mc() const { return mM; };
    std::size_t nnz() const { return mNNZ; };
    const mat_size_t &cn() const { return mCN; };
    const mat_size_t &vn() const { return mVN; };
    const vec_size_t &r() const { return mEdgeCN; };
    const vec_size_t &c() const { return mEdgeVN; };
    std::size_t nct() const { return mNCT; };
    std::size_t kct() const { return mKCT; };
    std::size_t mct() const { return mMCT; };
    const vec_size_t &puncture() const { return mPuncture; };
    const vec_size_t &shorten() const { return mShorten; };
    std::size_t max_dc() const { return mMaxDC; };

private:
    std::size_t mN;
    std::size_t mK;
    std::size_t mM;
    std::size_t mNNZ;
    mat_size_t mCN;       /* denotes the check neighbors, i.e. connected VN, for each check node as index in c/r; dimensions cn[mc][cw[i]] */
    mat_size_t mVN;       /* denotes the var neighbors, i.e., connected CN, for each variable node as index in c/r; dimensions vn[nc][vw[i]] */
    vec_size_t mEdgeCN;        /* non zero row indices; length nnz */
    vec_size_t mEdgeVN;        /* non zero check indices; length nnz */
    vec_size_t mPuncture; /* array pf punctured bit indices */
    vec_size_t mShorten;  /* array of shortened bit indices */
    std::size_t mNCT;     /* number of transmitted code bits */
    std::size_t mKCT;     /* number of transmitted information bits */
    std::size_t mMCT;     /* number of transmitted parity check bits */
    std::size_t mMaxDC;
};




} // namespace pgd