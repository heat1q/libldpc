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
        ldpc_decoder(const std::shared_ptr<ldpc_code> &code, const decoder_param &decoderParam);

        unsigned decode();
        bool is_codeword();

        //getter functions
        const vec_double_t &lv2c() const { return mLv2c; }
        const vec_double_t &lc2v() const { return mLc2v; }
        const vec_double_t &llr_in() const { return mLLRIn; }
        const vec_double_t &llr_out() const { return mLLROut; }

        const vec_bits_t &syndrome() const { return mSynd; }
        const vec_bits_t &estm_cw() const { return mCO; }

        static inline constexpr int sign(const double x);
        static inline constexpr double jacobian(const double x, const double y);
        static inline constexpr double minsum(const double x, const double y);

    private:
        std::shared_ptr<ldpc_code> mLdpcCode;

        vec_double_t mLv2c;
        vec_double_t mLc2v;

        // auxillary vectors for efficient CN update
        vec_double_t mExMsgF;
        vec_double_t mExMsgB;

        vec_double_t mLLRIn;
        vec_double_t mLLROut;

        vec_bits_t mSynd;
        vec_bits_t mCO;

        std::function<double(double, double)> mCNApprox;
        const decoder_param &mDecoderParam;
    };
} // namespace ldpc
