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
        friend class channel_bec;

    public:
        ldpc_decoder() = default;
        ldpc_decoder(const std::shared_ptr<ldpc_code> &code,
                     const decoder_param &decoderParam);
        virtual ~ldpc_decoder() = default;

        void set_approximation();

        int decode();
        bool is_codeword();

        // Set the input LLR
        void llr_in(const vec_double_t &llrIn) { mLLRIn = llrIn; }

        // Set the parameters
        void set_param(const decoder_param &decoderParam)
        {
            mDecoderParam = decoderParam;
            set_approximation();
        }

        //getter functions
        const decoder_param &param() const { return mDecoderParam; }

        const vec_double_t &lv2c() const { return mLv2c; }
        const vec_double_t &lc2v() const { return mLc2v; }
        const vec_double_t &llr_out() const { return mLLROut; }

        // The current estimated codeword
        const vec_bits_t &estimate() const { return mCO; }

        static inline constexpr int sign(const double x);
        static inline constexpr double jacobian(const double x, const double y);
        static inline constexpr double minsum(const double x, const double y);

    protected:
        std::shared_ptr<ldpc_code> mLdpcCode;

        vec_double_t mLv2c;
        vec_double_t mLc2v;

        // auxillary vectors for efficient CN update
        vec_double_t mExMsgF;
        vec_double_t mExMsgB;

        vec_double_t mLLRIn;
        vec_double_t mLLROut;

        vec_bits_t mCO;

        std::function<double(double, double)> mCNApprox;
        decoder_param mDecoderParam;
    };

    class ldpc_decoder_bec : public ldpc_decoder
    {
        friend class channel_bec;
    public:
        ldpc_decoder_bec() = default;
        ldpc_decoder_bec(const std::shared_ptr<ldpc_code> &code,
                         const decoder_param &decoderParam);
        virtual ~ldpc_decoder_bec() = default;

        int decode();

    private:
        std::vector<u8> mLv2c; 
        std::vector<u8> mLc2v;

        std::vector<u8> mLLRIn;
        std::vector<u8> mLLROut;
    };
} // namespace ldpc
