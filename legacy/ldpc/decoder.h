#pragma once

#include "ldpc.h"

namespace pgd
{
class ldpc_sim;

/**
 * @brief LDPC Decoder class
 * 
 */
class ldpc_decoder
{
public:
    ldpc_decoder(ldpc_code *pCode, ldpc_sim *pSim, std::int16_t pI, bool pEarlyTerm);

    void calc_llrs(const vec_double_t &y, double sigma2);

    std::int16_t decode();
    bool is_codeword_legacy();

    //getter functions
    std::size_t max_iter() const { return mMaxIter; }
    bool early_termination() const { return mEarlyTerm; }

    ldpc_code *ldpc() const { return mLdpcCode; }
    const vec_double_t &lv2c() const { return mLv2c; }
    const vec_double_t &lc2v() const { return mLc2v; }
    const vec_double_t &llr_in() const { return mLLRIn; }
    const vec_double_t &llr_out() const { return mLLROut; }

    const vec_bits_t &syndrome() const { return mSynd; }
    const vec_bits_t &estm_cw() const { return mCO; }

private:
    ldpc_code *mLdpcCode;
    ldpc_sim *mSim;

    vec_double_t mLv2c;
    vec_double_t mLc2v;
    vec_double_t mExMsgCN;

    vec_double_t mLLRIn;
    vec_double_t mLLROut;

    vec_bits_t mSynd;
    vec_bits_t mCO;

    std::int16_t mMaxIter;
    bool mEarlyTerm;
};
} // namespace pgd