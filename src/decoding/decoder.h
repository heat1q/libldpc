#pragma once

#include "../core/ldpc.h"

namespace ldpc
{
    constexpr int sign(const double x)
    {
        return (1 - 2 * static_cast<int>(std::signbit(x)));
    }

    constexpr double jacobian(const double x, const double y)
    {
        return sign(x) * sign(y) * std::min(std::abs(x), std::abs(y)) + std::log((1 + std::exp(-std::abs(x + y))) / (1 + std::exp(-std::abs(x - y))));
    }

    constexpr double minsum(const double x, const double y)
    {
        return sign(x) * sign(y) * std::min(std::abs(x), std::abs(y));
    }

    /**
    * @brief LDPC Decoder base class
    * 
    */
    template <typename T>
    class ldpc_decoder_base
    {
    public:
        ldpc_decoder_base() = default;
        ldpc_decoder_base(const std::shared_ptr<ldpc_code> &code,
                          const decoder_param &decoderParam)
            : mLdpcCode(code),
              mCNApprox(ldpc::jacobian),
              mCO(code->nc()),
              mLv2c(code->nnz()), mLc2v(code->nnz()),
              mExMsgF(code->max_dc()), mExMsgB(code->max_dc()),
              mLLRIn(code->nc()), mLLROut(code->nc())
        {
            set_param(decoderParam);
        }
        virtual ~ldpc_decoder_base() = default;

        virtual int decode();

        // Verifies whether mCO is a codeword or not
        bool is_codeword()
        {
            //calc syndrome
            bits_t s;
            for (int i = 0; i < mLdpcCode->mc(); i++)
            {
                s = 0;
                for (const auto &hj : mLdpcCode->H().row_neighbor()[i])
                {
                    s += mCO[hj.nodeIndex];
                }
                if (s != 0)
                {
                    return false;
                }
            }
            return true;
        }

        // Set the input LLR
        void set_llr_in(const std::vector<T> &in) { mLLRIn = in; }

        // Get the output LLR
        const std::vector<T> &llr_out() const { return mLLROut; }

        // Set the decoder parameters & update the CN approximation operation
        void set_param(const decoder_param &param)
        {
            mDecoderParam = param;
            if (mDecoderParam.type == std::string("BP_MS"))
            {
                mCNApprox = ldpc::minsum;
            }
        }

        // The current estimated codeword
        const vec_bits_t &estimate() const { return mCO; }

    protected:
        std::shared_ptr<ldpc_code> mLdpcCode;

        decoder_param mDecoderParam;

        // CN approximation operation
        std::function<T(T, T)> mCNApprox;

        // Estimated codeword
        vec_bits_t mCO;

        // auxillary vectors for efficient CN update
        std::vector<T> mLv2c;
        std::vector<T> mLc2v;

        std::vector<T> mExMsgF;
        std::vector<T> mExMsgB;

        std::vector<T> mLLRIn;
        std::vector<T> mLLROut;
    };

    /**
     * @brief Standard LDPC BP decoder
     * 
     */
    class ldpc_decoder : public ldpc_decoder_base<double>
    {
    public:
        friend class channel;
        friend class channel_awgn;
        friend class channel_bsc;

        ldpc_decoder() = default;
        ldpc_decoder(const std::shared_ptr<ldpc_code> &code,
                     const decoder_param &decoderParam);
        virtual ~ldpc_decoder() = default;

        int decode() override;
    };

    /**
     * @brief Simplified BP decoder for the BEC
     * 
     */
    class ldpc_decoder_bec : public ldpc_decoder_base<u8>
    {
    public:
        friend class channel_bec;

        ldpc_decoder_bec() = default;
        ldpc_decoder_bec(const std::shared_ptr<ldpc_code> &code,
                     const decoder_param &decoderParam);
        virtual ~ldpc_decoder_bec() = default;

        int decode() override;
    };
} // namespace ldpc
