#pragma once

#include "../core/ldpc.h"

namespace ldpc
{
    /**
    * @brief LDPC Decoder class
    * 
    */
    class ldpc_decoder
    {
        friend class channel;
        friend class channel_awgn;
        friend class channel_bsc;

    public:
        ldpc_decoder() = default;
        ldpc_decoder(const ldpc_code *code, const unsigned iter, const bool earlyTerm);

        unsigned decode();
        bool is_codeword();

        //getter functions
        const u64 max_iter() const { return mMaxIter; }
        const bool early_termination() const { return mEarlyTerm; }

        const ldpc_code *ldpc() const { return mLdpcCode; }
        const vec_double_t &lv2c() const { return mLv2c; }
        const vec_double_t &lc2v() const { return mLc2v; }
        const vec_double_t &llr_in() const { return mLLRIn; }
        const vec_double_t &llr_out() const { return mLLROut; }

        const vec_bits_t &syndrome() const { return mSynd; }
        const vec_bits_t &estm_cw() const { return mCO; }

    private:
        const ldpc_code *mLdpcCode;

        vec_double_t mLv2c;
        vec_double_t mLc2v;
        vec_double_t mExMsgCN;

        vec_double_t mLLRIn;
        vec_double_t mLLROut;

        vec_bits_t mSynd;
        vec_bits_t mCO;

        const unsigned mMaxIter;
        const bool mEarlyTerm;
    };
} // namespace ldpc
